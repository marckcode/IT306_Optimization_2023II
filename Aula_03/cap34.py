from pyomo.environ import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

m = AbstractModel()

# declaração de conjuntos e parâmetros
m.T = Set()                                 # conjunto de intervalos de tempo
m.PD = Param(m.T, domain=PositiveReals)     # demanda de potência ativa da residencia
m.PGD = Param(m.T, domain=NonNegativeReals) # geração distribuida da residencia
m.cOS = Param(m.T, domain=PositiveReals)    # custo pela compra de energia
m.PMXAE = Param(domain=NonNegativeReals)    # potência ativa máxima de injeção e extração do armazenador 
m.EMXAE = Param(domain=NonNegativeReals)    # capacidade máxima de armazenamento de energia 
m.EAE0 = Param(domain=NonNegativeReals)     # energia armazenada inicial
m.alpha = Param(domain=PositiveReals)       # eficiência do armazenador 
m.beta = Param(domain=PositiveReals)        # taxa de auto descarga do armazenador de energia
m.Pe = Param(domain=NonNegativeReals)       # potencia ativa da carga gerenciável. 
m.Te = Param(domain=NonNegativeReals)       # tempo de operação da carga gerenciável a ser ligada.
m.delta = Param(domain=PositiveReals)           # intervalo de tempo em horas    

# declaração das variáveis 
m.PS = Var(m.T, domain=NonNegativeReals)    # potência ativa fornecida pela rede de distribuição
m.PV = Var(m.T, domain=NonNegativeReals)    # potência ativa vendida ou entregue para a rede de distribuição
m.PAEi = Var(m.T, domain=NonNegativeReals)  # potência ativa de injeção do armazenador
m.PAEe = Var(m.T, domain=NonNegativeReals)  # potência ativa de extração do armazenador
m.EAE = Var(m.T, domain=Reals)              # energia armazenada 
m.eAE = Var(m.T, domain=Binary)             # estado de operação do armazenador   
m.we = Var(m.T, domain=Binary)              # horário de operação da carga gerenciável (1 operando, 0 cc)

# definição da função objetivo
# minimizar a energia total de consumo
def obj_func(m):
    return m.delta * sum(m.cOS[t] * m.PS[t] for t in m.T)
m.custo_total = Objective(rule=obj_func, sense=minimize)

# balanço de potência ativa na residencia
def balanço_potencia_ativa(m, t):
    return m.PS[t] + m.PGD[t] - m.PV[t] == m.PD[t] + m.PAEe[t] - m.PAEi[t] + m.Pe * m.we[t]
m.balanço_potencia_ativa = Constraint(m.T, rule=balanço_potencia_ativa)

def energia_armazenada(m, t):
    if t > 1:
        return m.EAE[t] == m.EAE[t-1] + m.alpha*m.delta*m.PAEe[t] \
                         - m.delta*m.PAEi[t]/m.alpha - m.beta*m.delta*m.EAE[t]
    else:
        return Constraint.Skip
m.energia_armazenada = Constraint(m.T, rule=energia_armazenada)

def energia_armazenada_inicial(m, t):
    if t == 1:
        return m.EAE[t] == m.EAE0 + m.alpha*m.delta*m.PAEe[t] \
                         - m.delta*m.PAEi[t]/m. alpha - m.beta*m.delta * m.EAE[t]
    else:
        return Constraint.Skip
m.energia_armazenada_inicial = Constraint(m.T, rule=energia_armazenada_inicial)

def operacao_carga_gerenciavel(m):
    return sum(m.delta * m.we[t] for t in m.T) == m.Te
m.operacao_carga_gerenciavel = Constraint(rule=operacao_carga_gerenciavel)

def tempo_operacao_carga_gerenciavel(m, t):
    if 2 <= t <= len(m.T) - m.Te/m.delta + 1:
        return sum(m.delta * m.we[k] for k in range(t, t + int(m.Te/m.delta))) >= \
                   m.Te * (m.we[t] - m.we[t - 1])
    else:
        return Constraint.Skip
m.tempo_operacao_carga_gerenciavel = Constraint(m.T, rule=tempo_operacao_carga_gerenciavel)

def tempo_operacao_carga_gerenciável_inicial(m, t):
    if t == 1:
        return sum(m.delta * m.we[k] for k in range(t, t + int(m.Te/m.delta))) >= m.Te * m.we[t]
    else:
        return Constraint.Skip
m.tempo_operacao_carga_gerenciável_inicial = Constraint(m.T, rule=tempo_operacao_carga_gerenciável_inicial)

def potência_ativa_máxima_injecao_armazenador(m, t):
    return m.PAEi[t] <= m.PMXAE * m.eAE[t]
m.potência_ativa_máxima_injecao_armazenador = Constraint(m.T, rule=potência_ativa_máxima_injecao_armazenador)

def potência_ativa_máxima_extracao_armazenador(m, t):
    return m.PAEe[t] <= m.PMXAE * (1 - m.eAE[t])
m.potência_ativa_máxima_extracao_armazenador = Constraint(m.T, rule=potência_ativa_máxima_extracao_armazenador)

def capacidade_máxima_armazenamento(m, t):
    return inequality(0, m.EAE[t], m.EMXAE)
m.capacidade_máxima_armazenamento = Constraint(m.T, rule=capacidade_máxima_armazenamento)

file = 'data/cap34.dat'
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
eAE= [instance.eAE[i].value for i in instance.eAE]
we = [instance.we[i].value for i in instance.we]

dt = [instance.delta*t for t in instance.T]
paee_paei = [instance.PAEe[t].value - instance.PAEi[t].value for t in instance.T]
pe_we = [instance.Pe*instance.we[t].value for t in instance.T]

display_df = {'PS': pot, 'PD': dem, 'PV': pvv, 'PGD': pgd, 
        'PAEi': paei, 'PAEe': paee, 'EAE': eae, 'eAE': eAE, 'we': we}
display_df = pd.DataFrame(display_df).map(format_value)

print(display_df)
print('\n')

# mostrar os resultados
print("Matriz[*]") 
matriz = []
for t in instance.T:
    delta_t = instance.delta * t
    row = [delta_t, instance.PS[t].value, instance.PD[t], instance.PGD[t], instance.PV[t].value,
           instance.PAEe[t].value-instance.PAEi[t].value, instance.EAE[t].value, 
           abs(instance.Pe*instance.we[t].value)]
    matriz.append(row)
   
columns = ['Delta*T', 'PS', 'PD', 'PGD', 'PV','PAEei', 'EAE', 'pe*we']
matriz = pd.DataFrame(matriz, columns=columns)
print(matriz)
print('\n')

t   = np.array(matriz.iloc[:, 0])
PS  = np.array(matriz.iloc[:, 1])
PD  = np.array(matriz.iloc[:, 2])
PGE = np.array(matriz.iloc[:, 3])
PV  = np.array(matriz.iloc[:, 4])
PAE = np.array(matriz.iloc[:, 5])
PE  = np.array(matriz.iloc[:, 7])


fig, ax = plt.subplots(figsize=(15, 8))

stacked_values = np.vstack([PGE - PV, PS])
bar_width = 0.8

ax.bar(t, stacked_values[0], width=bar_width, label='PGD-PV')
ax.bar(t, stacked_values[1], width=bar_width, bottom=stacked_values[0], label='PS')
ax.bar(t, -PV, color='green', label='PV')

ax.plot(t, PD, 'ko-', linewidth=2, markersize=5, label='PD')
ind = np.where(PE > 0)
ax.plot(t[ind], PE[ind]+PD[ind], 'ks-', linewidth=2, markersize=5, markerfacecolor='white', label='PE')

for i in range(len(t)):
    if PAE[i] < 0:
        ax.plot(t[i], PS[i] + PGE[i] - PV[i] - PAE[i] / 2, 'k^', linewidth=1, markersize=5, markerfacecolor='k')
    elif PAE[i] > 0:
        ax.plot(t[i], PS[i] + PGE[i] - PV[i] - PAE[i] / 2, 'kv', linewidth=1, markersize=5, markerfacecolor='k')
        

for i in range(len(t)):
    if PAE[i] < 0:
        ax.plot([t[i], t[i]], [PS[i] + PGE[i] - PV[i], PS[i] + PGE[i] - PV[i] - PAE[i]], 'k--', linewidth=1)
    elif PAE[i] > 0:
        ax.plot([t[i], t[i]], [PS[i] + PGE[i] - PV[i], PS[i] + PGE[i] - PV[i] - PAE[i]], 'k--', linewidth=1)

ax.legend(loc="upper left")
ax.set_xlabel('Tempo [h]', fontsize=12)
ax.set_xticks(np.arange(1, 25))
ax.set_ylabel('Potencia Ativa [kW]', fontsize=12)
fig.savefig('../imgs/cap34.png')

plt.show()
