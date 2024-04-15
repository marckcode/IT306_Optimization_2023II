from pyomo.environ import *
import pandas as pd

m = AbstractModel()

# declaração de conjuntos e parâmetros
m.SID = Set()           # Conjunto de siderúrgicas
m.FAB = Set()           # Conjunto de fábricas

m.COST = Param(m.SID, m.FAB, domain=NonNegativeReals)  # Custos de transporte por tonelada
m.PROD = Param(m.SID, domain=NonNegativeReals)         # Toneladas produzidas pelas siderúrgicas
m.DEMA = Param(m.FAB, domain=NonNegativeReals)         # Toneladas demandadas pelas fábricas

# Definir variáveis
m.TRANS = Var(m.SID, m.FAB, domain=NonNegativeReals)  # Toneladas a serem transportadas

# definição da função objetivo
# custo total de transporte
def func_obj(m):
    return sum(m.COST[s, f] * m.TRANS[s, f] for s in m.SID for f in m.FAB)
m.custo_total = Objective(rule=func_obj, sense=minimize)

# Definir restrições (por exemplo, demanda das fábricas)
def siderurgicas(m, s):
    return sum(m.TRANS[s, f] for f in m.FAB) == m.PROD[s]
m.siderurgicas = Constraint(m.SID, rule=siderurgicas)

def fabricas(m, f):
    return sum(m.TRANS[s, f] for s in m.SID) == m.DEMA[f]
m.fabricas = Constraint(m.FAB, rule=fabricas)

file = 'data/transporte.dat'
instance = m.create_instance(file)

path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
solver = SolverFactory("cplex", executable=path)  
results = solver.solve(instance, tee=False) # Mudar: tee=True -> show solver log

# transformar os dados a dataframe
# TRANS: 3x7 | SID: 3cols | FAB: 7rows
TRANS = [[instance.TRANS[s, f].value for s in instance.SID] for f in instance.FAB]
print('\n')
print('TRANS [*,*] (tr)')
df = pd.DataFrame(TRANS, index=instance.FAB, columns=instance.SID).sort_index()
print(df)
print('\n')
print(f"Custo total de transporte: {instance.custo_total()}")
print('\n')