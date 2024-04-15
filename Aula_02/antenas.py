from pyomo.environ import *
import pandas as pd

m = AbstractModel()

# declaração de conjuntos e parâmetros
m.MED = Set()  # Conjunto de medidores
m.ANT = Set()  # Conjunto de antenas
m.BAI = Set()  # Conjunto de bairros

m.DIST = Param(m.MED, m.ANT, domain=NonNegativeReals)  # Matriz de distâncias
m.DMAX = Param(domain=NonNegativeReals)                # Distância máxima de leitura
m.CIO  = Param(m.ANT, domain=NonNegativeReals)         # Custo de instalação e operação
m.BB   = Param(m.ANT, domain=NonNegativeReals)         # Bairro

# declaração das variáveis
m.LOCA = Var(m.ANT, domain=Binary)  # localizacao de antenas
m.COVE = Var(m.MED, m.ANT, domain=Binary)  # matriz de covertura  

# definição da função objetivo
# custo total de transporte
def obj_func(m):
    return sum(m.CIO[a] * m.LOCA[a] for a in m.ANT)
m.numero_antenas = Objective(rule=obj_func, sense=minimize)

# definição das restrições
def aloca_antena(m, med, ant):
    return m.COVE[med, ant] <= m.LOCA[ant]
m.aloca_antena = Constraint(m.MED, m.ANT, rule=aloca_antena)

def distancia_maxima(m, med, ant):
    return m.DIST[med, ant] * m.COVE[med, ant] <= m.DMAX
m.distancia_maxima = Constraint(m.MED, m.ANT, rule=distancia_maxima)

def covertura_minima(m, med):
    return sum(m.COVE[med, ant] for ant in m.ANT) >= 2
m.covertura_minima = Constraint(m.MED, rule=covertura_minima)

def posicao_bairro(m, b):
    return sum(m.LOCA[ant] for ant in m.ANT if m.BB[ant] == b) >= 1
m.posicao_bairro = Constraint(m.BAI, rule=posicao_bairro)


file = 'data/antenas.dat'
instance = m.create_instance(file)
path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
solver = SolverFactory('cplex', executable=path)
results = solver.solve(instance, tee=True) # Mudar: tee=False -> skip solver log

print(f"\nObjective: {instance.numero_antenas()}") # print objective

# mudar para o dataframe para melhor visualização dos dados
loc_values = [instance.LOCA[l].value for l in instance.LOCA]
df_loc = pd.DataFrame(loc_values, index=instance.ANT, columns=['LOC'])
df_loc = df_loc.astype(int)
print('\n')
print('LOCA [*] := ')
print(df_loc)

# TRANS: 8x10 | ANT: 8 | MED: 10
ant_values = [[instance.COVE[f, s].value for s in instance.ANT] for f in instance.MED]
df_ant = pd.DataFrame(ant_values, index=instance.MED, columns=instance.ANT)
df_ant = df_ant.astype(int)
print('\n')
print('COVE [*,*]')
print(df_ant)
print('\n')