from pyomo.environ import *

model = AbstractModel()

# declaração de conjuntos e parâmetros
model.P = Set(dimen=1)          # conjunto de pontos tabelados
model.x = Param(model.P)        # coordenada x tabelados
model.y = Param(model.P)        # coordenada y tabelados 

# declaração das variáveis
model.ycalc = Var(model.P)      # coordenada y calculada
model.m = Var()                 # inclinação da reta ou coeficiente linear
model.b = Var()                 # coeficiente constante 

# definição da função objetivo
# minimizar a soma dos quadrados das diferenças entre os valores 
# tabelados e os valores obtidos pela aproximação
def func_obj(model):
    return sum((model.ycalc[p] - model.y[p])**2 for p in model.P)
model.minimos_quadrados = Objective(rule=func_obj, sense=minimize)

def calculo_ycalc(model, p):
    return model.ycalc[p] == model.m * model.x[p] + model.b
model.calculo_ycalc = Constraint(model.P, rule=calculo_ycalc)

file = 'data/min_quadrados.dat' 
instance = model.create_instance(file)

path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
solver = SolverFactory("cplex", executable=path)
results = solver.solve(instance, tee=False) # Mudar: tee=True -> show solver log

print("Status: " + str(results.solver.status))
print("Result: " + str(results.solver.termination_condition))
print('\n')
for i in instance.P:
    print('Var ycalc[%d] = %2.4f' % (i, value(instance.ycalc[i])))
print('Var m: %2.4f' % value(instance.m))
print('Var b: %2.4f' % value(instance.b))
print('\n')
print('Objective: %2.4f' % instance.minimos_quadrados())
print('\n')