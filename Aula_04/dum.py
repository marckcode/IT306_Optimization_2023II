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
m.cIPA = Param()        # custo da capacidade da potência ativa máxima do armazenador [R$/kW].
m.cIEA = Param()        # custo da capacidade máxima de armazenamento de energia [R$/kWh].
m.cOS = Param(m.T)      # custo pela compra de energia da rede de distribuição [R$/kWh]   
m.cOT = Param()         # custo de operação do sistema de geração convencional [R$/kWh]
m.cCC = Param()         # custo do corte ou racionamento de carga [R$/kWh]
m.EAE0 = Param()        # energia armazenada inicial [kWh] 
m.alpha = Param()       # eficiência do armazenador de energia [%] 
m.beta = Param()        # taxa de auto descarga do armazenador de energia [%]
m.delta = Param()       # intervalo de tempo em horas [h]
m.D = Param()           # tempo de duração da microrrede operando de forma isolada em horas [h]
m.p = Param(m.C)        # probabilidade do cenário de operação [%]

# declaração das variáveis
m.PGDmax = Var(domain=NonNegativeReals, 
               bounds=(0, m.IPGDmax))           # tamanho da geração distribuída renovável da microrrede [kW] 
m.PTmax = Var(domain=NonNegativeReals, 
               bounds=(0, m.IPTmax))            # tamanho da geração convencional da microrrede [kW]
m.PAEmax = Var(domain=NonNegativeReals)         # potência ativa máxima de injeção e extração do armazenador [kW]
m.EAEmax = Var(domain=NonNegativeReals)         # capacidade máxima de armazenamento de energia [kWh]
m.PS = Var(m.T, m.C, domain=NonNegativeReals)   # potência ativa fornecida pela rede de distribuição [kW]
m.PT = Var(m.T, m.C, domain=NonNegativeReals)   # potência ativa fornecida pela geração convencional [kW]
m.xD = Var(m.T, m.C, domain=NonNegativeReals, 
                     bounds=(0, 1))             # porcentagem de corte ou racionamento de carga [%] 
m.PAEi = Var(m.T, m.C, domain=NonNegativeReals) # potência ativa de injeção do armazenador [kW]
m.PAEe = Var(m.T, m.C, domain=NonNegativeReals) # potência ativa de extração do armazenador [kW]
m.EAE = Var(m.T, m.C, domain=NonNegativeReals)  # energia armazenada [kWh]

# definição da função objetivo
# minimizar o custo do dimensionamento
def custo(m):
    return m.cIGD*m.PGDmax + m.cIT * m.PTmax + m.cIPA * m.PAEmax + m.cIEA * m.EAEmax + \
                    365 * sum(m.p[c] * m.delta * m.cOS[t] * m.PS[t, c] for t in m.T for c in m.C) + \
                    365 * sum(m.p[c] * m.delta * m.cCC * m.PD[t] * m.xD[t, c] for t in m.T for c in m.C) + \
                    365 * sum(m.p[c] * m.delta * m.cOT * m.PT[t, c] for t in m.T for c in m.C)
m.custo = Objective(rule=custo, sense=minimize)

# balanço de potência ativa 
def balanco_potencia_ativa(m, t, c):
    return m.PS[t, c] + m.PT[t, c] + m.fGD[t] * m.PGDmax == m.PD[t] * (1 - m.xD[t, c]) + m.PAEe[t, c] - m.PAEi[t, c]
m.balanco_potencia_ativa = Constraint(m.T, m.C, rule=balanco_potencia_ativa)

def capacidade_da_subestacao(m, t, c):
    return m.PS[t, c] <= m.PSmax
m.capacidade_da_subestacao = Constraint(m.T, m.C, rule=capacidade_da_subestacao)

def capacidade_do_gerador_convencional(m, t, c):
    return m.PT[t, c] <= m.PTmax
m.capacidade_do_gerador_convencional = Constraint(m.T, m.C, rule=capacidade_do_gerador_convencional)

def potência_ativa_máxima_injecao_armazenador(m, t, c):
    return m.PAEi[t, c] <= m.PAEmax
m.potência_ativa_máxima_injecao_armazenador = Constraint(m.T, m.C, rule=potência_ativa_máxima_injecao_armazenador)

def potência_ativa_máxima_extracao_armazenador(m, t, c):
    return m.PAEe[t, c] <= m.PAEmax
m.potência_ativa_máxima_extracao_armazenador = Constraint(m.T, m.C, rule=potência_ativa_máxima_extracao_armazenador)

def energia_armazenada(m, t, c):
    if t > 1:
        return m.EAE[t, c] == m.EAE[t-1, c] + m.alpha * m.delta * m.PAEe[t, c] - \
                              m.delta * m.PAEi[t, c] / m.alpha - m.beta * m.delta * m.EAE[t, c]
    else:
        return Constraint.Skip
m.energia_armazenada = Constraint(m.T, m.C, rule=energia_armazenada)

def energia_armazenada_inicial(m, t, c):
    if t == 1:
        return m.EAE[t, c] == m.EAE0 + m.alpha * m.delta * m.PAEe[t, c] - \
                              m.delta * m.PAEi[t, c] / m.alpha - m.beta * m.delta * m.EAE[t, c]
    else:
        return Constraint.Skip
m.energia_armazenada_inicial = Constraint(m.T, m.C, rule=energia_armazenada_inicial)

def capacidade_máxima_armazenamento(m, t, c):
    return m.EAE[t, c] <= m.EAEmax
m.capacidade_máxima_armazenamento = Constraint(m.T, m.C, rule=capacidade_máxima_armazenamento)

def operacao_da_subestacao_sob_contingencia(m, t, c):
    if c != 0 and c <= t <= min(len(m.T), c + m.D/m.delta - 1):
        return m.PS[t, c] == 0
    else:
        return Constraint.Skip
m.operacao_da_subestacao_sob_contingencia = Constraint(m.T, m.C, rule=operacao_da_subestacao_sob_contingencia)


file = 'data/dum.dat'
instance = m.create_instance(file)

path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
solver = SolverFactory('cplex', executable=path)
results = solver.solve(instance, tee=False)
print('\nCusto total: %.4f' % instance.custo())

# imprime os resultados
print(f"\nPGDmax = {instance.PGDmax.value}")
print(f"\ncIGD * PGDmax = {int(instance.cIGD.value * instance.PGDmax.value)}")

print(f"\nPTmax = {instance.PTmax.value}")
print(f"\ncIT * PTmax = {int(instance.cIT.value * instance.PTmax.value)}")

print(f"\nPAEmax = {instance.PAEmax.value}")
print(f"\ncIPA * PAEmax = {int(instance.cIPA.value * instance.PAEmax.value)}")

print(f"\nEAEmax = {instance.EAEmax.value}")
print(f"\ncIEA * EAEmax = {int(instance.cIEA.value * instance.EAEmax.value)}")

res1 = 0
for t in instance.T:
    for c in instance.C:
        res1 += instance.p[c] * instance.delta * instance.cOS[t] * instance.PS[t, c].value
print(f"\n365*(sum'{'t in T'}' sum'{'c in C'}' p[c]*delta*cOS[t]*PS[t,c]) = {ceil(365 * res1)}") 

res2 = 0
for t in instance.T:
    for c in instance.C:
        res2 += instance.p[c] * instance.delta * instance.cOT * instance.PT[t, c].value
print(f"\n365 * sum'{'t in T'}' sum'{'c in C'}' (p[c]*delta*cOT*PT[t,c]) = {round(365 * res2, 1)}") 

res3 = 0
for t in instance.T:
    for c in instance.C:
        res3 += instance.p[c] * instance.delta * instance.cCC * instance.PD[t] *instance.xD[t, c].value
print(f"\n365*(sum'{'t in T'}' sum'{'c in C'}' p[c]*delta*cCC*PD[t]*xD[t,c])) = {int(365 * res3)}") 

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

print("\nPS [*,*]")
ps_values = {c: [instance.PS[t, c].value for t in instance.T] for c in instance.C}
ps_df = pd.DataFrame(ps_values).map(format_value)
print(ps_df)
print('\n')

print("xD [*,*]")
xd_values = {c: [instance.xD[t, c].value for t in instance.T] for c in instance.C}
xd_df = pd.DataFrame(xd_values).map(format_value)
print(xd_df)
print('\n')

print("PT [*,*]")
pt_values = {c: [instance.PT[t, c].value for t in instance.T] for c in instance.C}
pt_df = pd.DataFrame(pt_values).map(format_value)
print(pt_df)
print('\n')

print("PAEi [*,*]")
paei_values = {c: [instance.PAEi[t, c].value for t in instance.T] for c in instance.C}
paei_df = pd.DataFrame(paei_values).map(format_value)
print(paei_df)
print('\n')

print("PAEe [*,*]")
paee_values = {c: [instance.PAEe[t, c].value for t in instance.T] for c in instance.C}
paee_df = pd.DataFrame(paee_values).map(format_value)
print(paee_df)
print('\n')

print("EAE [*,*]")
eae_values = {c: [instance.EAE[t, c].value for t in instance.T] for c in instance.C}
eae_df = pd.DataFrame(eae_values).map(format_value)
print(eae_df)
print('\n')

