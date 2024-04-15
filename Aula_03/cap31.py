from pyomo.environ import *
import pandas as pd

m = AbstractModel()

# declaração de conjuntos e parâmetros
m.T = Set()                                     # conjunto de intervalos de tempo
m.PD = Param(m.T, domain=PositiveReals)         # demanda de potência ativa da residencia
m.cOS = Param(m.T, domain=PositiveReals)        # custo pela compra de energia
m.delta = Param(domain=PositiveReals)           # intervalo de tempo em horas    

# declaração das variáveis 
m.PS = Var(m.T, domain=NonNegativeReals)         # potência ativa fornecida pela rede de distribuição

# definição da função objetivo
# minimizar a energia total de consumo
def obj_func(m):
    return m.delta * sum(m.cOS[t] * m.PS[t] for t in m.T)
m.custo_total = Objective(rule=obj_func, sense=minimize)

# balanço de potência ativa na residencia
def balanço_potencia_ativa(m, t):
    return m.PS[t] == m.PD[t]
m.balanço_potencia_ativa = Constraint(m.T, rule=balanço_potencia_ativa)

file = 'data/cap31.dat'
instance = m.create_instance(file)

path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
solver = SolverFactory('cplex', executable=path)
results = solver.solve(instance, tee=False)

print('\nObjective: %.4f' % instance.custo_total())

# Definir uma função para customizar o dataframe
def format_value(value, decimals=4):
    if abs(value) == 0 or abs(value) == 1:
        return "{:.0f}".format(abs(value))
    elif int(abs(value)) == value:
        return str(int(value))
    elif str(value).count('.') == 1 and len(str(value).split('.')[1]) == 1:
        return "{:.1f}".format(value)  # Display with one decimal place
    else:
        return "{:.{}f}".format(value, decimals)

# create the data frame
pot = [instance.PS[i].value for i in instance.PS]
dem = [instance.PD[i] for i in instance.PD]
data = {'PS': pot, 'PD': dem}
df = pd.DataFrame(data).map(format_value)

print(df)
print('\n')











