from pyomo.environ import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

m = AbstractModel()

# usar servidor NEOS para acelerar os resultados 
os.environ['NEOS_EMAIL'] = 'm272292@dac.unicamp.br'

# Definição dos conjuntos e parâmetros
m.S = Set()         # Conjunto de estados
m.T = Set()         # Conjunto de discretizações

m.D = Param(m.S)        # ID do dispositivo
m.Pdisp = Param(m.S, domain=NonNegativeReals)    # Potência Ativa dos estados
m.Ptotal = Param(m.T, domain=NonNegativeReals)   # Potência total da leitura do medidor
m.previo = Param(m.S)                             # estado previo
m.MD = Param(m.S)                                 # tempo minimo de operação

# Definição das variáveis
m.x = Var(m.S, m.T, domain=Binary)          # estado do dispositivo ‘e’ para o instante de tempo ‘t’
m.delta_P = Var(m.T)                        # erro para o instante de tempo ‘t’
m.up = Var(m.S, m.T, domain=Binary)         # variáveis auxiliares 
m.dw = Var(m.S, m.T, domain=Binary)         # variáveis auxiliares 

# Função objetivo
def erro(m):
    return sum(m.delta_P[t] for t in m.T)
m.erro = Objective(rule=erro, sense=minimize)

# Definição das restrições
def diferenca_combinatoria_1(m, t):
    return m.Ptotal[t] - sum(m.Pdisp[s] * m.x[s, t] for s in m.S) <= m.delta_P[t]
m.diferenca_combinatoria_1 = Constraint(m.T, rule=diferenca_combinatoria_1)

def diferenca_combinatoria_2(m, t):
    return m.Ptotal[t] - sum(m.Pdisp[s] * m.x[s, t] for s in m.S) >= -m.delta_P[t]
m.diferenca_combinatoria_2 = Constraint(m.T, rule=diferenca_combinatoria_2)

def sobreposicao_estados(m, t, d):
    return sum(m.x[s, t] for s in m.S if m.D[s] == d) <= 1
m.sobreposicao_estados = Constraint(m.T, range(1, 4), rule=sobreposicao_estados)

def sequencia_estados_1(m, s, t):
    if t > 1:
        return m.x[s, t] - m.x[s, t-1] == m.up[s, t] - m.dw[s, t]
    else:
        return Constraint.Skip
m.sequencia_estados_1 = Constraint(m.S, m.T, rule=sequencia_estados_1)

def sequencia_estados_2(m, s, t):
    return m.up[s, t] + m.dw[s, t] <= 1
m.sequencia_estados_2 = Constraint(m.S, m.T, rule=sequencia_estados_2)

def sequencia_estados_3(m, s, t):
    if m.previo[s] > 0:
        return m.up[s, t] == m.dw[m.previo[s], t]
    else:
        return Constraint.Skip
m.sequencia_estados_3 = Constraint(m.S, m.T, rule=sequencia_estados_3)

def tempo_minimo_operacao_rule(m, s, t):
    if t in range(2, len(m.T) - m.MD[s] + 2):
        return sum(m.x[s, k] for k in range(t, t + m.MD[s])) >= m.MD[s] * (m.x[s, t] - m.x[s, t-1])
    else:
        return Constraint.Skip
m.tempo_minimo_operacao = Constraint(m.S, m.T, rule=tempo_minimo_operacao_rule)


data = DataPortal() # conjunto de dados, cuando tem-se muitos archivos .dat
data.load(filename="data/NILM.dat", model=m)
data.load(filename="data/dados_medidor.dat", model=m)
instance = m.create_instance(data)

# path = "/mnt/c/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
# solver = SolverFactory('cplex', executable=path)
# results = solver.solve(instance, tee=False)

solver_manager = SolverManagerFactory('neos')
results = solver_manager.solve(instance, solver="cplex")
print('\nErro: %.2f' % instance.erro())

delta_P = []
for t in instance.T:
    delta_P.append([int(t), instance.delta_P[t].value])
delta_P = pd.DataFrame(delta_P, columns=['Index', 'deltaP'])

x = {}
for t in instance.T:
    x[t] = [int(instance.x[s, t].value) for s in instance.S]
x = pd.DataFrame(x.values(), index=x.keys(), columns=['E1', 'E2', 'E3', 'E4', 'E5', 'E6']) 

Ptotal = []
for t in instance.T:
    Ptotal.append([int(t), instance.Ptotal[t]])
nPtotalorg = pd.DataFrame(Ptotal, columns=['Index', 'Ptotal'])

Pdisp = []
for s in instance.S:
    Pdisp.append([int(s), instance.Pdisp[s]])
Pdisp = pd.DataFrame(Pdisp, columns=['Index', 'Pdisp'])

# ================================ IMPRECÁO DE DADOS
# Ler dados originais
M = pd.read_csv('data/ground_truth.csv', names=["Val", "refr", "wash", "stove"])

refr = M.iloc[:, 0]
wash = M.iloc[:, 1]
stove = M.iloc[:, 2]
all = M.iloc[:, 3]
time = list(range(1, len(refr)+1))

# Ler arquivos com dados
DELTA = delta_P
X = x
ESTADO = Pdisp

# Vetores
TS = list(X.index.values)
Ptotal = list(np.dot(X.iloc[:, :], ESTADO.iloc[:, 1]).squeeze())

# Matriz para associacáo de estados aos respectivos dispositivos
states = [1, 1, 2, 2, 2, 3]

num_rows = max(states)
num_cols = len(states)
F = np.zeros((num_rows, num_cols))

for i in range(num_cols):        
    for j in range(num_rows):    
        if j + 1 == states[i]:
            F[j, i] = 1
            
K = np.dot(F, np.diag(ESTADO.iloc[:, 1]))
Pdisp = np.dot(X.iloc[:, :], K.T).squeeze()


# Gráfico dos resultados
fig, ax = plt.subplots(2, 1, figsize=(10,6))

time = np.array(time)
refr = np.array(refr)
wash = np.array(wash)
stove = np.array(stove)

# Dados Ground Truth
bw = 1
ax[0].bar(time,  refr, width=bw, label='Geladeira')
ax[0].bar(time,  wash, width=bw, bottom=refr, label='Máquina de lavar')
ax[0].bar(time, stove, width=bw, bottom=refr + wash, label='Forno elétrico')

ax[0].set_xlim(0, 700)
ax[0].set_ylim(0, 3500)
ax[0].legend()
ax[0].set_xlabel('t (min)')
ax[0].set_ylabel('Potência Ativa [W]')
ax[0].set_title('Dados Originais')

# Os resultados
ax[1].bar(time, Pdisp[:, 0], width=bw)
ax[1].bar(time, Pdisp[:, 1], width=bw, bottom=Pdisp[:, 0])
ax[1].bar(time, Pdisp[:, 2], width=bw, bottom=Pdisp[:, 0] + Pdisp[:, 1])

ax[1].set_xlim(0, 700)
ax[1].set_ylim(0, 3500)
ax[1].set_xlabel('Amostra (delta = 15s)')
ax[1].set_ylabel('Potência Ativa [W]')
ax[1].set_title('Modelo Matemático')

plt.tight_layout()
plt.show()

fig.savefig('../imgs/NILM.png')

