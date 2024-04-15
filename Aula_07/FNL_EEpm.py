from pyomo.environ import *
import pandas as pd


m = AbstractModel()

# declaração de conjuntos e parâmetros
m.ob   = Set()                      # conjunto de nós
m.ol   = Set(dimen=2)               # conjunto de circuitos
m.obV  = Set()                      # conjunto de medidas de magnitude de tensão
m.olP  = Set(dimen=2)               # conjunto de medidas de fluxo de potencia ativa
m.obPS = Set()                      # conjunto de medidas de potencia ativa da subestação

m.obout = Set(m.ob, within=m.ob, initialize=[])
m.obin  = Set(m.ob, within=m.ob, initialize=[])

def Populate_In_and_Out(m):
    for i, j in m.ol:
        m.obin[j].add(i)
        m.obout[i].add(j)
m.In_n_Out = BuildAction(rule=Populate_In_and_Out)

m.Tb     = Param(m.ob)                       # tipo de barra 0: carga, 1: SE
m.PDprev = Param(m.ob)                       # potência ativa de demanda prevista no nó i
m.QDprev = Param(m.ob)                       # potência reativa de demanda prevista no nó i
m.R      = Param(m.ol, mutable=True)         # resistência no circuito ij
m.X      = Param(m.ol, mutable=True)         # reatância no circuito ij
m.Z2     = Param(m.ol, mutable=True)         # impedância no circuito ij ao quadrado
m.Imax   = Param(m.ol)                       # limite máximo da magnitude de corrente no circuito ij
m.Vnom   = Param(mutable=True)               # magnitude da tensão nominal  
m.Vmax   = Param(mutable=True)               # magnitude de tensão amáxima
m.Vmin   = Param(mutable=True)               # magnitude de tensão mínima
m.Vmed   = Param(m.obV)                      # medidas de magnitude de tensão
m.Pmed   = Param(m.olP)                      # medidas de fluxo de potencia ativa
m.PSmed  = Param(m.obPS)                     # medidas de potencia ativa da subestação
m.WV     = Param(m.obV)                      # ponderações das medidas de magnitude de tensão
m.WP     = Param(m.olP)                      # ponderações das medidas de fluxo de potencia ativa
m.WPS    = Param(m.obPS)                     # ponderações das medidas de potencia ativa da subestação

# declaração das variáveis
m.V     = Var(m.ob)                          # variável que representa a magnitude de tensão
m.PS    = Var(m.ob)                          # potência ativa fornecida pela subestação no nó i
m.QS    = Var(m.ob)                          # potência reativa fornecida pela subestação no nó i
m.I     = Var(m.ol)                          # variável que representa o fluxo de corrente
m.P     = Var(m.ol)                          # fluxo de potência ativa no circuito ij
m.Q     = Var(m.ol)                          # fluxo de potência reativa no circuito ij
m.eV    = Var(m.obV)                         # resíduos da magnitude de tensão   
m.eP    = Var(m.olP)                         # resíduos do fluxo de potencia ativa
m.ePS   = Var(m.obPS)                        # resíduos da potencia ativa da subestação
m.alpha = Var(initialize=1.0)                # ajuste da previsão da demanda às medidas existentes 

# definição da função objetivo
def minimos_quadrados(m):
    return sum(m.eV[k]*m.WV[k]*m.eV[k] for k in m.obV) + \
           sum(m.eP[k,l]*m.WP[k,l]*m.eP[k,l] for k, l in m.olP) + \
           sum(m.ePS[k]*m.WPS[k]*m.ePS[k] for k in m.obPS) 
m.minimos_quadrados = Objective(rule=minimos_quadrados, sense=minimize)

# definição das restrições
def calculo_residuos_tensão(m, k):
    return m.Vmed[k] == m.V[k] + m.eV[k]
m.calculo_residuos_tensão = Constraint(m.obV, rule=calculo_residuos_tensão)

def calculo_residuos_potencia_ativa(m, k, l):
    return m.Pmed[k, l] == m.P[k, l] + m.eP[k, l]
m.calculo_residuos_potencia_ativa = Constraint(m.olP, rule=calculo_residuos_potencia_ativa)

def calculo_residuos_injecao_ativa(m, k):
    return m.PSmed[k] == m.PS[k] + m.ePS[k]
m.calculo_residuos_injecao_ativa = Constraint(m.obPS, rule=calculo_residuos_injecao_ativa)

def balance_potencia_activa(m, i):
    return sum(m.P[j, i] for j in m.obin[i]) - \
           sum(m.P[i, j] + m.R[i, j]*m.I[i, j]**2 for j in m.obout[i]) + m.PS[i] == m.alpha * m.PDprev[i]
m.balance_potencia_activa = Constraint(m.ob, rule=balance_potencia_activa)

def balance_potencia_reativa(m, i):
    return sum(m.Q[j, i] for j in m.obin[i]) - \
           sum(m.Q[i, j] + m.X[i, j]*m.I[i, j]**2 for j in m.obout[i]) + m.QS[i] == m.alpha * m.QDprev[i]
m.balance_potencia_reativa = Constraint(m.ob, rule=balance_potencia_reativa)

def queda_magnitude_tensao(m, i, j):
    return m.V[i]**2 - 2*(m.R[i, j]*m.P[i, j] + m.X[i, j]*m.Q[i, j]) - m.Z2[i,j] * m.I[i, j]**2 - m.V[j]**2 == 0
m.queda_magnitude_tensao = Constraint(m.ol, rule=queda_magnitude_tensao)

def calculo_magnitude_corrente(m, i, j):
    return m.V[j]**2 * m.I[i,j]**2 == m.P[i,j]**2 + m.Q[i,j]**2 
m.calculo_magnitude_corrente = Constraint(m.ol, rule=calculo_magnitude_corrente)

def limite_magnitude_corrente(m, i, j):
    return inequality(0, m.I[i, j], m.Imax[i, j])
m.limite_magnitude_corrente = Constraint(m.ol, rule=limite_magnitude_corrente)

def limite_magnitude_tensão(m, i):
    return inequality(m.Vmin, m.V[i], m.Vmax)
m.limite_magnitude_tensão = Constraint(m.ob, rule=limite_magnitude_tensão)

file = "data/sistema34nospm.dat"
inst = m.create_instance(file)

inst.Vnom = 11.0/sqrt(3)   
inst.Vmax = 11.0/sqrt(3)            
inst.Vmin = 9.0/sqrt(3)

# fixar a injeção de potencia ativa nos nós de carga iguais a zero
for i in inst.ob:
    if inst.Tb[i] == 0:
        inst.PS[i].fix(0)
        
# fixar a injeção de potencia reativa nos nós de carga iguais a zero
for i in inst.ob:
    if inst.Tb[i] == 0:
        inst.QS[i].fix(0)
        
# nomalização e calculo da impedancia ao quadrado
for i, j in inst.ol:
    inst.R[i, j]  = inst.R[i, j] / 1000      # de ohm para kohm
    inst.X[i, j]  = inst.X[i, j] / 1000      # de ohm para kohm
    inst.Z2[i, j] = inst.R[i, j]**2 + inst.X[i, j]**2
    
# ponto inicial
for i in inst.ob:
    inst.V[i]  = inst.Vnom
    inst.PS[i] = 0
    inst.QS[i] = 0
    
for i, j in inst.ol:
    inst.I[i, j] = 0
    inst.P[i, j] = 0
    inst.Q[i, j] = 0

for i in inst.obV:
    inst.eV[i] = 0
    inst.V[i]  = inst.Vmed[i]

for i, j in inst.olP:
    inst.eP[i, j] = 0
    inst.P[i, j]  = inst.Pmed[i, j]

    
for i in inst.obPS:
    inst.ePS[i] = 0
    inst.PS[i]  = inst.PSmed[i]

solver = SolverFactory('ipopt')
results = solver.solve(inst, tee=True)
print('Minimos Quadrados: %.4f' % inst.minimos_quadrados())


# definir uma funcão de customização
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
print()
df1 = []
for i in inst.ob:
    row = [i, inst.V[i].value, inst.PS[i].value, inst.QS[i].value]
    df1.append(row)
        
columns = ['i', 'V', 'PS', 'QS']
df1 = pd.DataFrame(df1, columns=columns).map(format_value)
print(df1)

print()
df2 = []
for i, j in inst.ol:
    row = [i, j, inst.I[i,j].value, inst.P[i,j].value, inst.Q[i,j].value]
    df2.append(row)
        
columns = ['i', 'j', 'I', 'P', 'Q']
df2 = pd.DataFrame(df2, columns=columns).map(format_value)
print(df2)

# display eV, eP, ePS, alpha;
print()
print("alpha = ", inst.alpha.value)
print()
df3 = []
for i in inst.obV:
    row = [i, inst.eV[i].value]
    df3.append(row)
        
columns = ['i', 'eV']
df3 = pd.DataFrame(df3, columns=columns).map(format_value)
print(df3)

print()
df4 = []
for i in inst.obPS:
    row = [i, inst.ePS[i].value]
    df4.append(row)
        
columns = ['i', 'ePS']
df4 = pd.DataFrame(df4, columns=columns).map(format_value)
print(df4)

print()