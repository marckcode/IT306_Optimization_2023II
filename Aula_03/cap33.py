from pyomo.environ import *
import pandas as pd

m = AbstractModel()

# declaração de conjuntos e parâmetros
m.T = Set()                                     # conjunto de intervalos de tempo
m.PD = Param(m.T, domain=PositiveReals)         # demanda de potência ativa da residencia
m.PGD = Param(m.T, domain=NonNegativeReals)     # geração distribuida da residencia
m.cOS = Param(m.T, domain=PositiveReals)        # custo pela compra de energia
m.PMXAE = Param(domain=NonNegativeReals)        # potência ativa máxima de injeção e extração do armazenador 
m.EMXAE = Param(domain=NonNegativeReals)        # capacidade máxima de armazenamento de energia 
m.EAE0 = Param(domain=NonNegativeReals)         # energia armazenada inicial
m.alpha = Param(domain=PositiveReals)           # eficiência do armazenador 
m.beta = Param(domain=NonNegativeReals)         # taxa de auto descarga do armazenador de energia
m.delta = Param(domain=PositiveReals)           # intervalo de tempo em horas    

# declaração das variáveis 
m.PS = Var(m.T, domain=NonNegativeReals)        # potência ativa fornecida pela rede de distribuição
m.PV = Var(m.T, domain=NonNegativeReals)        # potência ativa vendida ou entregue para a rede de distribuição
m.PAEi = Var(m.T, domain=NonNegativeReals)      # potência ativa de injeção do armazenador
m.PAEe = Var(m.T, domain=NonNegativeReals)      # potência ativa de extração do armazenador
m.EAE = Var(m.T, domain=Reals)                  # energia armazenada
m.eAE = Var(m.T, domain=Boolean)                # estado de operação do armazenador   

# definição da função objetivo
# minimizar a energia total de consumo
def obj_func(m):
    return m.delta * sum(m.cOS[t] * m.PS[t] for t in m.T)
m.custo_total = Objective(rule=obj_func, sense=minimize)

# balanço de potência ativa na residencia
def balanço_potencia_ativa(m, t):
    return m.PS[t] + m.PGD[t] - m.PV[t] == m.PD[t] + m.PAEe[t] - m.PAEi[t]
m.balanço_potencia_ativa = Constraint(m.T, rule=balanço_potencia_ativa)

def energia_armazenada(m, t):
    if t > 1:
        return m.EAE[t] == m.EAE[t-1] + m.alpha*m.delta*m.PAEe[t] - m.delta*m.PAEi[t]/m.alpha - m.beta*m.delta*m.EAE[t]
    else:
        return Constraint.Skip
m.energia_armazenada = Constraint(m.T, rule=energia_armazenada)

def energia_armazenada_inicial(m, t):
    if t == 1:
        return m.EAE[t] == m.EAE0 + m.alpha*m.delta*m.PAEe[t] - m.delta*m.PAEi[t]/m. alpha - m.beta*m.delta * m.EAE[t]
    else:
        return Constraint.Skip
m.energia_armazenada_inicial = Constraint(m.T, rule=energia_armazenada_inicial)

def potência_ativa_máxima_injecao_armazenador(m, t):
    return m.PAEi[t] <= m.PMXAE * m.eAE[t]
m.potência_ativa_máxima_injecao_armazenador = Constraint(m.T, rule=potência_ativa_máxima_injecao_armazenador)

def potência_ativa_máxima_extracao_armazenador(m, t):
    return m.PAEe[t] <= m.PMXAE * (1 - m.eAE[t])
m.potência_ativa_máxima_extracao_armazenador = Constraint(m.T, rule=potência_ativa_máxima_extracao_armazenador)

def capacidade_máxima_armazenamento(m, t):
    return inequality(0, m.EAE[t], m.EMXAE)
m.capacidade_máxima_armazenamento = Constraint(m.T, rule=capacidade_máxima_armazenamento)

file = 'data/cap33.dat'
instance = m.create_instance(file)

path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
solver = SolverFactory('cplex', executable=path)
results = solver.solve(instance, tee=False)
print('\nCusto total: %.4f' % instance.custo_total())

# Definir uma função para customizar o dataframe
def format_value(value, decimals=4):
    if abs(value) == 0 or abs(value) == 1:
        return "{:.0f}".format(abs(value))
    elif int(abs(value)) == value:
        return str(int(value))
    elif str(value).count('.') == 1 and len(str(value).split('.')[1]) == 1:
        return "{:.1f}".format(value) 
    else:
        return "{:.{}f}".format(value, decimals)

# criar o dataframe
pot = [instance.PS[i].value for i in instance.PS]
dem = [instance.PD[i] for i in instance.PD]
pvv = [instance.PV[i].value for i in instance.PV]
pgd = [instance.PGD[i] for i in instance.PGD]
paei = [instance.PAEi[i].value for i in instance.PAEi]
paee = [instance.PAEe[i].value for i in instance.PAEe]
eae = [instance.EAE[i].value for i in instance.EAE]
Eae = [instance.eAE[i].value for i in instance.eAE]

data = {'PS': pot, 'PD': dem, 'PV': pvv, 'PGD': pgd, 'PAEi': paei, 'PAEe': paee, 'EAE': eae, 'eAE': Eae}
df = pd.DataFrame(data).map(format_value)

print(df)
print('\n')