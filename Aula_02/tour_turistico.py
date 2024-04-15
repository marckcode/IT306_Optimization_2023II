from pyomo.environ import *
import pandas as pd

# Create a abstract model
m = AbstractModel()

# declaração de conjuntos e parâmetros
m.STUR = Set()                          # conjunto de sitios turisticos

m.SATIS = Param(m.STUR, domain=PositiveIntegers)  # Satisfação esperada
m.CENTR = Param(m.STUR, domain=NonNegativeReals)  # Custo da entrada
m.TVISI = Param(m.STUR, domain=NonNegativeReals)  # Tempo de visita
m.ORCAM = Param(domain=NonNegativeReals)          # Orçamento
m.TEMAX = Param(domain=NonNegativeReals)          # Tempo máximo de visita

# declaração das variáveis
m.DESCI = Var(m.STUR, domain=Binary)                                   # decisao entrar ou nao

# definição da função objetivo
# custo total de transporte
def func_obj(m):
    return sum(m.SATIS[s]*m.DESCI[s] for s in m.STUR)
m.satisfacao_esperada = Objective(rule=func_obj, sense=maximize)

def orcamento(m):
    return sum(m.CENTR[s]*m.DESCI[s] for s in m.STUR) <= m.ORCAM
m.orcamento = Constraint(rule=orcamento)

def tempo_maximo(m):
    return sum(m.TVISI[s]*m.DESCI[s] for s in m.STUR) <= m.TEMAX
m.tempo_maximo = Constraint(rule=tempo_maximo)


file = 'data/tour_turistico.dat'
instance = m.create_instance(file)

path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
solver = SolverFactory("cplex", executable=path)
solver.options['mipgap'] = 0
results = solver.solve(instance, tee=False) # Mudar: tee=True -> show solver log

DESCI = {}
for s in instance.STUR:
    i = int(instance.DESCI[s].value)
    DESCI[s] = i
    
# transformar os dados a dataframe
print('\n')
df = pd.DataFrame(list(DESCI.items()), columns=['City_Tour', 'DESCI'])
print('DESCI [*] :=')
print(df)

print('\n')
print('Gasto:     ', sum(instance.CENTR[s]*instance.DESCI[s].value for s in instance.STUR))
print('T. Visita: ', sum(instance.TVISI[s]*instance.DESCI[s].value for s in instance.STUR))
print(f"Satisfacao: {instance.satisfacao_esperada()}")
print('\n')






