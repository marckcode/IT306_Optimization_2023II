from pyomo.environ import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

m = AbstractModel()

# declaração de conjuntos e parâmetros
m.T = Set()                                       # conjunto de intervalos de tempo
m.AE = Set()                                      # conjunto de armazenadores de energia
m.CG = Set()                                      # conjunto de cargas gerenciáveis
# declaração de params    
m.PD = Param(m.T, within=PositiveReals)           # demanda de potência ativa da residencia
m.PGD = Param(m.T, within=NonNegativeReals)       # geração distribuida da residencia
m.cOS = Param(m.T, within=PositiveReals)          # custo pela compra de energia
m.PMXAE = Param(m.AE)                             # potência ativa máxima de injeção e extração do armazenador 
m.EMXAE = Param(m.AE)                             # capacidade máxima de armazenamento de energia 
m.EAE0 = Param(m.AE)                              # energia armazenada inicial
m.alpha = Param(m.AE)                             # eficiência do armazenador 
m.delta = Param(within=PositiveReals)             # intervalo de tempo em horas    
m.beta = Param(m.AE)                              # taxa de auto descarga do armazenador de energia
m.Pe = Param(m.CG)                                # potencia ativa da carga gerenciável 
m.Te = Param(m.CG)                                # tempo de operação da carga gerenciável a ser ligada no dia 
m.Ti = Param(m.CG)                                # tempo de operação inicial da carga gerenciável
m.Tf = Param(m.CG)                                # tempo de operação final da carga gerenciável
# declaração das variáveis 
m.PS = Var(m.T, within=NonNegativeReals)          # potência ativa fornecida pela rede de distribuição
m.PV = Var(m.T, within=NonNegativeReals)          # potência ativa vendida ou entregue para a rede de distribuição
m.PAEi = Var(m.T, m.AE, within=NonNegativeReals)  # potência ativa de injeção do armazenador
m.PAEe = Var(m.T, m.AE, within=NonNegativeReals)  # potência ativa de extração do armazenador
m.EAE = Var(m.T, m.AE)                            # energia armazenada 
m.eAE = Var(m.T, m.AE, domain=Binary)             # estado de operação do armazenador   
m.we = Var(m.T, m.CG, domain=Binary)              # horário de operação da carga gerenciável (1 operando, 0 cc)

# definição da função objetivo
# minimizar a energia total de consumo
def obj_func(m):
    return m.delta * sum(m.cOS[t] * m.PS[t] for t in m.T)
m.custo_total = Objective(rule=obj_func, sense=minimize)

# balanço de potência ativa na residencia
def balanço_potencia_ativa(m, t):
    return m.PS[t] + m.PGD[t] - m.PV[t] == m.PD[t] \
           + sum(m.PAEe[t, a] for a in m.AE) - sum(m.PAEi[t, a] for a in m.AE) \
           + sum(m.Pe[c] * m.we[t, c] for c in m.CG)
m.balanço_potencia_ativa = Constraint(m.T, rule=balanço_potencia_ativa)

def balanco_energetico(m, t, a):
    if t > 1:
        return m.EAE[t, a] == m.EAE[t-1, a] + m.alpha[a]*m.delta*m.PAEe[t, a] \
             - m.delta*m.PAEi[t, a] / m.alpha[a] - m.beta[a]*m.delta*m.EAE[t, a]
    else:
        return Constraint.Skip
m.balanco_energetico = Constraint(m.T, m.AE, rule=balanco_energetico)

def balanco_energetico_zero(m, t, a):
    if t == 1:
        return m.EAE[t, a] == m.EAE0[a] + m.alpha[a]*m.delta*m.PAEe[t, a] \
             - m.delta*m.PAEi[t, a] / m. alpha[a] \
             - m.beta[a]*m.delta * m.EAE[t, a]
    else:
        return Constraint.Skip
m.balanco_energetico_zero = Constraint(m.T, m.AE, rule=balanco_energetico_zero)

def carga_gerenciavel_ligado(m, c):
    return sum(m.delta * m.we[t, c] for t in m.T) == m.Te[c]
m.carga_gerenciavel_ligado = Constraint(m.CG, rule=carga_gerenciavel_ligado)

def carga_gerenciavel_duracao(m, t, c):
    if 2 <= t <= len(m.T) - m.Te[c]/m.delta + 1:
        return sum(m.delta * m.we[k, c] for k in range(t, t + int(m.Te[c]/m.delta))) \
                >= m.Te[c] * (m.we[t, c] - m.we[t-1, c])
    else:
        return Constraint.Skip
m.carga_gerenciavel_duracao = Constraint(m.T, m.CG, rule=carga_gerenciavel_duracao)

def carga_gerenciavel_duracao_inicial(m, c, t):
    if t == 1:
        return sum(m.delta * m.we[k, c] for k in range(t, t + int(m.Te[c]/m.delta))) >= \
                   m.Te[c] * m.we[t, c]
    else:
        return Constraint.Skip
m.carga_gerenciavel_duracao_inicial = Constraint(m.CG, m.T, rule=carga_gerenciavel_duracao_inicial)

def carga_gerenciavel_desligado(m, t, c):
    if t < m.Ti[c] / m.delta or t > m.Tf[c] / m.delta:
        return m.we[t, c] == 0
    else:
        return Constraint.Skip
m.carga_gerenciavel_desligado = Constraint(m.T, m.CG, rule=carga_gerenciavel_desligado)

def potência_ativa_máxima_injecao_armazenador(m, t, a):
    return m.PAEi[t,a] <= m.PMXAE[a] * m.eAE[t,a]
m.potência_ativa_máxima_injecao_armazenador = Constraint(m.T, m.AE, rule=potência_ativa_máxima_injecao_armazenador)

def potência_ativa_máxima_extracao_armazenador(m, t, a):
    return m.PAEe[t,a] <= m.PMXAE[a] * (1 - m.eAE[t,a])
m.potência_ativa_máxima_extracao_armazenador = Constraint(m.T, m.AE, rule=potência_ativa_máxima_extracao_armazenador)

def capacidade_energia(m, t, a):
    return inequality(0, m.EAE[t, a], m.EMXAE[a])
m.capacidade_energia = Constraint(m.T, m.AE, rule=capacidade_energia)

file = 'data/geer.dat'
instance = m.create_instance(file)

path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
solver = SolverFactory('cplex', executable=path)
results = solver.solve(instance, tee=False)
print('Custo total: %.4f' % instance.custo_total())

# salvar os dados nos dataframes 
matriz = []
for t in instance.T:
    delta_t = instance.delta * t
    row = [delta_t, instance.PS[t].value, instance.PD[t], instance.PGD[t], instance.PV[t].value]
    matriz.append(row)

columns = ['Delta*T', 'PS', 'PD', 'PGD', 'PV']
matriz = pd.DataFrame(matriz, columns=columns) # df1

matrizPAE = []
for t in instance.T:
    row = [instance.PAEe[t, a].value - instance.PAEi[t, a].value for a in instance.AE]
    matrizPAE.append(row)
        
columns = [f'PAE_{a}' for a in instance.AE]
matrizPAE = pd.DataFrame(matrizPAE, columns=columns) # df2

matrizEAE = []
for t in instance.T:
    row = [instance.EAE[t, a].value for a in instance.AE]
    matrizEAE.append(row)
        
columns = [f'EAE_{a}' for a in instance.AE]
matrizEAE = pd.DataFrame(matrizEAE, columns=columns) # df3


matrizPWE = []
for t in instance.T:
    row = [abs(instance.Pe[c]*instance.we[t,c].value) for c in instance.CG]
    matrizPWE.append(row)
        
columns = [f'PWE_{c}' for c in instance.CG]
matrizPWE = pd.DataFrame(matrizPWE, columns=columns)


t   = np.array(matriz.iloc[:, 0])
PS  = np.array(matriz.iloc[:, 1])
PD  = np.array(matriz.iloc[:, 2])
PGE = np.array(matriz.iloc[:, 3])
PV  = np.array(matriz.iloc[:, 4])

PAE = (matrizPAE.T).sum(axis=0).T
PE = (matrizPWE.T).sum(axis=0).T

fig, ax = plt.subplots(figsize=(15, 8))

stacked_values = np.vstack([PGE - PV, PS])
bar_width = 0.2

ax.bar(t, stacked_values[0], width=bar_width, label='Geraçao Distr. Injet. a Residencia')
ax.bar(t, stacked_values[1], width=bar_width, bottom=stacked_values[0], label='Rede de Distribuçao')
ax.bar(t, -PV, color='green', label='Geraçao Distr. Vendida ou Entregue a Rede')
ax.plot(t, PD, 'ko-', linewidth=2, markersize=5, label='Demanda da Residencia')
ind = np.where(PE == 0)
PED = PE + PD
PED.iloc[ind] = np.nan
ax.plot(t, PED, 'ks-', linewidth=2, markersize=5, markerfacecolor='white', label='Carga Gerenciavel')

PAEaux = []
for i in range(len(t)):
    if PAE[i] < 0:
        PAEaux.append([t[i], PS[i] + PGE[i] - PV[i] - PAE[i] / 2])
PAEaux = np.array(PAEaux)

ax.plot(PAEaux[:, 0], PAEaux[:, 1], 'k^', linewidth=1, markersize=5, 
                                           markerfacecolor='k', label='Armazenador Descarregando')

PAEaux = []
for i in range(len(t)):
    if PAE[i] > 0:
        PAEaux.append([t[i], PS[i] + PGE[i] - PV[i] - PAE[i] / 2])
PAEaux = np.array(PAEaux)

ax.plot(PAEaux[:, 0], PAEaux[:, 1], 'kv', linewidth=1, markersize=5, 
                                           markerfacecolor='k', label='Armazenador Carregando')

for i in range(len(t)):
    if PAE[i] < 0:
        ax.plot([t[i], t[i]], [PS[i] + PGE[i] - PV[i], PS[i] + PGE[i] - PV[i] - PAE[i]], 
                                        'k--', linewidth=1)
    elif PAE[i] > 0:
        ax.plot([t[i], t[i]], [PS[i] + PGE[i] - PV[i], PS[i] + PGE[i] - PV[i] - PAE[i]], 
                                        'k--', linewidth=1)

ax.legend(loc="upper left")

ax.set_xlabel('Tempo [h]', fontsize=12)
ax.set_ylabel('Potencia Ativa [kW]', fontsize=12)

fig.savefig('../imgs/geer.png')
plt.show()







