from pyomo.environ import *
import pandas as pd

m = AbstractModel()

# declaraçao de conjuntos e parâmetros
m.T = Set()             # conjunto de intervalos de tempo
m.C = Set()             # conjunto de cenários de operação 

m.PD = Param(m.T)       # demanda de potência ativa da microrrede [kW]
m.fGD = Param(m.T)      # informação estatística diária normalizada da geração distribuída 
                        # renovável disponível em função a seu tamanho
m.PSmax = Param()       # capacidade da rede de distribuição que fornece energia à 
                        # microrrede em modo conectado [kW]
m.IPGDmax = Param()     # limite máximo para o tamanho da geração distribuída renovável [kW]
m.IPTmax = Param()      # limite máximo para o tamanho da geração convencional [kW]    
m.cIGD = Param()        # custo da capacidade da geração distribuída renovável [R$/kW]
m.cIT = Param()         # custo da capacidade da geração convencional [R$/kW] 
m.cOS = Param(m.T)      # custo pela compra de energia da rede de distribuição [R$/kWh]   
m.cOT = Param()         # custo de operação do sistema de geração convencional [R$/kWh]
m.cCC = Param()         # custo do corte ou racionamento de carga [R$/kWh]
m.delta = Param()       # intervalo de tempo em horas [h]    
m.D = Param()           # tempo de duração da microrrede operando de forma isolada em horas [h]
m.p = Param(m.C)        # probabilidade do cenário de operação [%]

# declaração das variáveis
m.PGDmax = Var(within=NonNegativeReals, 
               bounds=(0, m.IPGDmax))           # tamanho da geração distribuída renovável da microrrede [kW] 
m.PTmax = Var(within=NonNegativeReals, 
               bounds=(0, m.IPTmax))            # tamanho da geração convencional da microrrede [kW]
m.PS = Var(m.T, m.C, within=NonNegativeReals)   # potência ativa fornecida pela rede de distribuição [kW]
m.PT = Var(m.T, m.C, within=NonNegativeReals)   # potência ativa fornecida pela geração convencional [kW]
m.xD = Var(m.T, m.C, within=NonNegativeReals, 
                     bounds=(0, 1))             # porcentagem de corte ou racionamento de carga [%] 

# definição da função objetivo
# minimizar o custo do dimensionamento
def custo(m):
    return (m.cIGD*m.PGDmax +  m.cIT * m.PTmax + 
                365 * sum(m.p[c] * m.delta * m.cOS[t] * m.PS[t, c] for t in m.T for c in m.C) + 
                365 * sum(m.p[c] * m.delta * m.cOT * m.PT[t, c] for t in m.T for c in m.C) +
                365 * sum(m.p[c] * m.delta * m.cCC * m.PD[t] * m.xD[t, c] for t in m.T for c in m.C))
m.custo = Objective(rule=custo, sense=minimize)

# balanço de potência ativa 
def balanco_potencia_ativa(m, t, c):
    return m.PS[t, c] + m.PT[t, c] + m.fGD[t] * m.PGDmax == m.PD[t] * (1 - m.xD[t, c])
m.balanco_potencia_ativa = Constraint(m.T, m.C, rule=balanco_potencia_ativa)

def capacidade_da_subestacao(m, t, c):
    return m.PS[t, c] <= m.PSmax
m.capacidade_da_subestacao = Constraint(m.T, m.C, rule=capacidade_da_subestacao)

def capacidade_do_gerador_convencional(m, t, c):
    return m.PT[t, c] <= m.PTmax
m.capacidade_do_gerador_convencional = Constraint(m.T, m.C, rule=capacidade_do_gerador_convencional)

def operacao_da_subestacao_sob_contingencia(m, t, c):
    if c != 0 and c <= t <= min(len(m.T), c + m.D/m.delta - 1):
        return m.PS[t, c] == 0
    else:
        return Constraint.Skip
m.operacao_da_subestacao_sob_contingencia = Constraint(m.T, m.C, rule=operacao_da_subestacao_sob_contingencia)

file = 'data/cap42.dat'
instance = m.create_instance(file)

path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
solver = SolverFactory('cplex', executable=path)
results = solver.solve(instance, tee=False)
print('\nCusto total: %.4f' % instance.custo())

# Impressão dos resultados
print(f"PGDmax = {instance.PGDmax.value}")
print(f"cIGD * PGDmax = {int(instance.cIGD.value * instance.PGDmax.value)}")
print(f"PTmax = {instance.PTmax.value}")
print(f"cIT*PTmax = {int(instance.cIT.value * instance.PTmax.value)}")

res1 = 0
for t in instance.T:
    for c in instance.C:
        res1 += instance.p[c] * instance.delta * instance.cOS[t] * instance.PS[t, c].value
print(f"365*(sum'{'t in T'}' sum'{'c in C'}' p[c]*delta*cOS[t]*PS[t,c]) = {ceil(365 * res1)}") 

res2 = 0
for t in instance.T:
    for c in instance.C:
        res2 += instance.p[c] * instance.delta * instance.cOT * instance.PT[t, c].value
print(f"365 * sum'{'t in T'}' sum'{'c in C'}' (p[c]*delta*cOT*PT[t,c]) = {round(365 * res2, 1)}") 

res3 = 0
for t in instance.T:
    for c in instance.C:
        res3 += instance.p[c] * instance.delta * instance.cCC * instance.PD[t] *instance.xD[t, c].value
print(f"365*(sum'{'t in T'}' sum'{'c in C'}' p[c]*delta*cCC*PD[t]*xD[t,c])) = {int(365 * res3)}") 

print('\n')
# Definir uma função para customizar o dataframe
def format_value(value, decimals=4):
    if abs(value) == 0 or abs(value) == 1: 
        return "{:.0f}".format(value)
    elif int(abs(value)) == value:
        return str(int(value))
    elif str(value).count('.') == 1 and len(str(value).split('.')[1]) == 1:
        return "{:.1f}".format(value)
    else:
        return "{:.{}f}".format(value, decimals)

print('PS [*,*]')
ps_values = {c: [instance.PS[t, c].value for t in instance.T] for c in instance.C}
ps_df = pd.DataFrame(ps_values).map(format_value)
print(ps_df)
print('\n')

print('xD [*,*]')
xd_values = {c: [instance.xD[t, c].value for t in instance.T] for c in instance.C}
xd_df = pd.DataFrame(xd_values).map(format_value)
print(xd_df)
print('\n')

print('PT [*,*]')
pt_values = {c: [instance.PT[t, c].value for t in instance.T] for c in instance.C}
pt_df = pd.DataFrame(pt_values).map(format_value)
print(pt_df)
print('\n')
