from pyomo.environ import *

m = AbstractModel()

# declaração de conjuntos e parâmetros
m.PROD = Set()                                       # produtos
m.DISP = Param(domain=NonNegativeReals)          # horas disponíveis na semana
m.VELC = Param(m.PROD, domain=PositiveReals)     # toneladas de fabricação por hora
m.LUCR = Param(m.PROD)                               # lucro por toneladas
m.MAXT = Param(m.PROD, domain=NonNegativeReals)  # toneladas máximas de fabricação

ub = {'bandas':6000, 'bobinas':4000}

def fb(instance, i):
    return (0, ub[i])

# declaração das variáveis
m.FABR = Var(m.PROD, within=NonNegativeReals, bounds=fb)

# definição da função objetivo
# lucro total de produção
def func_obj(m):
    return sum(m.LUCR[p] * m.FABR[p] for p in m.PROD)
m.lucro_total = Objective(rule=func_obj, sense=maximize)

# definição das restrições
# o número total de horas utilizadas por todos os produtos 
# não pode exceder as horas disponíveis da fabrica
def tempo(m):
    return sum((1.0/m.VELC[p])*m.FABR[p] for p in m.PROD) <= m.DISP
m.tempo = Constraint(rule=tempo)

file = 'data/producao.dat'
instance = m.create_instance(file)

path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
solver = SolverFactory("cplex", executable=path)
results = solver.solve(instance, tee=False) # Mudar: tee=True -> show solver log

# Mostra os resultados
print('\n')
print('FABR [*] :=')
for p in instance.PROD:
    print(f"{p} = {instance.FABR[p].value}")
print('\n')
print(f"Lucro Total: {instance.lucro_total()}")
print('\n')








