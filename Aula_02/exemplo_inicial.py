from pyomo.environ import *

m = ConcreteModel()

# declaração das variáveis
m.x = Var()
m.y = Var()

# definição da função objetivo
def func_obj(m):
    return 25*m.x + 30*m.y
m.lucro_total = Objective(rule=func_obj, sense=maximize)

# definição das restrições
def tempo(m):
    return (1/200)*m.x + (1/140)*m.y <= 40
m.tempo = Constraint(rule=tempo)

def limite_bandas(m):
    return inequality(0, m.x, 6000)
m.limite_bandas = Constraint(rule=limite_bandas)

def limite_bobinas(m):
    return inequality(0, m.y, 4000)
m.limite_bobinas = Constraint(rule=limite_bobinas)

path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
solver = SolverFactory("cplex", executable=path)
results = solver.solve(m, tee=False)
print("Status: " + str(results.solver.status))
print("Result: " + str(results.solver.termination_condition))
print("\n")
print('Var x: ', value(m.x))
print('Var y: ', value(m.y))
print('Minimization OF: ' + str(m.lucro_total()))
print('\n')