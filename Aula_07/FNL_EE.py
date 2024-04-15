from pyomo.environ import *
import pandas as pd

m = AbstractModel()

# declaração de conjuntos e parâmetros
m.ob   = Set()                      # conjunto de nós
m.ol   = Set(dimen=2)               # conjunto de circuitos
m.obV  = Set()                      # conjunto de medidas de magnitude de tensão
m.olI  = Set(dimen=2)               # conjunto de medidas de magnitude de corrente
m.olP  = Set(dimen=2)               # conjunto de medidas de fluxo de potencia ativa
m.olQ  = Set(dimen=2)               # conjunto de medidas de fluxo de potencia reativa
m.obPS = Set()                      # conjunto de medidas de potencia ativa da subestação
m.obQS = Set()                      # conjunto de medidas de potencia reativa da subestação
m.obPD = Set()                      # conjunto de medidas de potencia ativa da demanda
m.obQD = Set()                      # conjunto de medidas de potencia reativa da demanda

m.obout = Set(m.ob, within=m.ob, initialize=[])
m.obin  = Set(m.ob, within=m.ob, initialize=[])

def Populate_In_and_Out(m):
    for i, j in m.ol:
        m.obin[j].add(i)
        m.obout[i].add(j)
m.In_n_Out = BuildAction(rule=Populate_In_and_Out)


m.Tb    = Param(m.ob)                       # tipo de barra 0: carga, 1: SE
m.R     = Param(m.ol, mutable=True)         # resistência no circuito ij
m.X     = Param(m.ol, mutable=True)         # reatância no circuito ij
m.Z2    = Param(m.ol, mutable=True)         # impedância no circuito ij ao quadrado
m.Imax  = Param(m.ol)                       # limite máximo da magnitude de corrente no circuito ij
m.Vnom  = Param(mutable=True)               # magnitude da tensão nominal  
m.Vmax  = Param(mutable=True)               # magnitude de tensão amáxima
m.Vmin  = Param(mutable=True)               # magnitude de tensão mínima
m.Vmed  = Param(m.obV)                      # medidas de magnitude de tensão
m.Imed  = Param(m.olI)                      # medidas de magnitude de corrente
m.Pmed  = Param(m.olP)                      # medidas de fluxo de potencia ativa
m.Qmed  = Param(m.olQ)                      # medidas de fluxo de potencia reativa
m.PSmed = Param(m.obPS)                     # medidas de potencia ativa da subestação
m.QSmed = Param(m.obQS)                     # medidas de potencia reativa da subestação
m.PDmed = Param(m.obPD)                     # medidas de potencia ativa da demanda
m.QDmed = Param(m.obQD)                     # medidas de potencia reativa da demanda
m.WV    = Param(m.obV)                      # ponderações das medidas de magnitude de tensão
m.WI    = Param(m.olI)                      # ponderações das medidas de magnitude de corrente
m.WP    = Param(m.olP)                      # ponderações das medidas de fluxo de potencia ativa
m.WQ    = Param(m.olQ)                      # ponderações das medidas de fluxo de potencia reativa
m.WPS   = Param(m.obPS)                     # ponderações das medidas de potencia ativa da subestação
m.WQS   = Param(m.obQS)                     # ponderações das medidas de potencia reativa da subestação
m.WPD   = Param(m.obPD)                     # ponderações das medidas de potencia ativa da demanda
m.WQD   = Param(m.obQD)                     # ponderações das medidas de potencia reativa da demanda


# declaração das variáveis
m.V   = Var(m.ob)                           # variável que representa a magnitude de tensão
m.PD  = Var(m.ob, domain=NonNegativeReals)  # potência ativa de demanda no nó i
m.QD  = Var(m.ob, domain=NonNegativeReals)  # potência reativa de demanda no nó i
m.PS  = Var(m.ob, domain=NonNegativeReals)  # potência ativa fornecida pela subestação no nó i
m.QS  = Var(m.ob, domain=NonNegativeReals)  # potência reativa fornecida pela subestação no nó i
m.I   = Var(m.ol)                           # variável que representa o fluxo de corrente
m.P   = Var(m.ol)                           # fluxo de potência ativa no circuito ij
m.Q   = Var(m.ol)                           # fluxo de potência reativa no circuito ij
m.eV  = Var(m.obV)                          # resíduos da magnitude de tensão 
m.eI  = Var(m.olI)                          # resíduos da magnitude de tensão     
m.eP  = Var(m.olP)                          # resíduos do fluxo de potencia ativa
m.eQ  = Var(m.olQ)                          # resíduos do fluxo de potencia reativa
m.ePS = Var(m.obPS)                         # resíduos da potencia ativa da subestação
m.eQS = Var(m.obQS)                         # resíduos da potencia reativa da subestação
m.ePD = Var(m.obPD)                         # resíduos da demanda de potencia ativa
m.eQD = Var(m.obQD)                         # resíduos da demanda de potencia reativa

# definição da função objetivo
def minimos_quadrados(m):
    return (sum(m.eV[k]*m.WV[k]*m.eV[k] for k in m.obV) + sum(m.eI[k,l]*m.WI[k,l]*m.eI[k,l] for k, l in m.olI) + 
            sum(m.eP[k,l]*m.WP[k,l]*m.eP[k,l] for k, l in m.olP) + sum(m.eQ[k,l]*m.WQ[k,l]*m.eQ[k,l] for k, l in m.olQ) + 
            sum(m.ePS[k]*m.WPS[k]*m.ePS[k] for k in m.obPS) + sum(m.eQS[k]*m.WQS[k]*m.eQS[k] for k in m.obQS) + 
            sum(m.ePD[k]*m.WPD[k]*m.ePD[k] for k in m.obPD) + sum(m.eQD[k]*m.WQD[k]*m.eQD[k] for k in m.obQD))
m.minimos_quadrados = Objective(rule=minimos_quadrados, sense=minimize)

# definição das restrições
def calculo_residuos_tensão(m, k):
    return m.Vmed[k] == m.V[k] + m.eV[k]
m.calculo_residuos_tensão = Constraint(m.obV, rule=calculo_residuos_tensão)

def calculo_residuos_corrente(m, k, l):
    return m.Imed[k, l] == m.I[k, l] + m.eI[k, l]
m.calculo_residuos_corrente = Constraint(m.olI, rule=calculo_residuos_corrente)

def calculo_residuos_potencia_ativa(m, k, l):
    return m.Pmed[k, l] == m.P[k, l] + m.eP[k, l]
m.calculo_residuos_potencia_ativa = Constraint(m.olP, rule=calculo_residuos_potencia_ativa)

def calculo_residuos_potencia_reativa(m, k, l):
    return m.Qmed[k, l] == m.Q[k, l] + m.eQ[k, l]
m.calculo_residuos_potencia_reativa = Constraint(m.olQ, rule=calculo_residuos_potencia_reativa)

def calculo_residuos_injecao_ativa(m, k):
    return m.PSmed[k] == m.PS[k] + m.ePS[k]
m.calculo_residuos_injecao_ativa = Constraint(m.obPS, rule=calculo_residuos_injecao_ativa)

def calculo_residuos_injecao_reativa(m, k):
    return m.QSmed[k] == m.QS[k] + m.eQS[k]
m.calculo_residuos_injecao_reativa = Constraint(m.obQS, rule=calculo_residuos_injecao_reativa)

def calculo_residuos_demanda_ativa(m, k):
    return m.PDmed[k] == m.PD[k] + m.ePD[k]
m.calculo_residuos_demanda_ativa = Constraint(m.obPD, rule=calculo_residuos_demanda_ativa)

def calculo_residuos_demanda_reativa(m, k):
    return m.QDmed[k] == m.QD[k] + m.eQD[k]
m.calculo_residuos_demanda_reativa = Constraint(m.obQD, rule=calculo_residuos_demanda_reativa)

def balance_potencia_activa(m, i):
    return (sum(m.P[j, i] for j in m.obin[i]) - 
           sum(m.P[i, j] + m.R[i, j]*m.I[i, j]**2 for j in m.obout[i]) + m.PS[i] == m.PD[i])
m.balance_potencia_activa = Constraint(m.ob, rule=balance_potencia_activa)

def balance_potencia_reativa(m, i):
    return (sum(m.Q[j, i] for j in m.obin[i]) - 
           sum(m.Q[i, j] + m.X[i, j]*m.I[i, j]**2 for j in m.obout[i]) + m.QS[i] == m.QD[i])
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



file = "data/sistema34nos.dat"
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
    inst.V[i] = inst.Vnom
    inst.PD[i] = 0
    inst.QD[i] = 0
    inst.PS[i] = 0
    inst.QS[i] = 0
    
for i, j in inst.ol:
    inst.I[i, j] = 0
    inst.P[i, j] = 0
    inst.Q[i, j] = 0

for i in inst.obV:
    inst.eV[i] = 0
    inst.V[i] = inst.Vmed[i]

for i, j in inst.olI:
    inst.eI[i, j] = 0
    inst.I[i, j] = inst.Imed[i, j]
    
for i, j in inst.olP:
    inst.eP[i, j] = 0
    inst.P[i, j] = inst.Pmed[i, j]
    
for i, j in inst.olQ:
    inst.eQ[i, j] = 0
    inst.Q[i, j] = inst.Qmed[i, j]
    
for i in inst.obPS:
    inst.ePS[i] = 0
    inst.PS[i] = inst.PSmed[i]
    
for i in inst.obQS:
    inst.eQS[i] = 0
    inst.QS[i] = inst.QSmed[i]
    
for i in inst.obPD:
    inst.ePD[i] = 0
    inst.PD[i] = inst.PDmed[i]
    
for i in inst.obQD:
    inst.eQD[i] = 0
    inst.QD[i] = inst.QDmed[i]


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
V  = [inst.V[i].value for i in inst.ob]
PD = [inst.PD[i].value for i in inst.ob]
QD = [inst.QD[i].value for i in inst.ob]
PS = [inst.PS[i].value for i in inst.ob]
QS = [inst.QS[i].value for i in inst.ob]
I  = [inst.I[i, j].value for i, j in inst.ol]
P  = [inst.P[i, j].value for i, j in inst.ol]
Q  = [inst.Q[i, j].value for i, j in inst.ol]
ind = [(i, j) for i, j in inst.ol]
eV  = [inst.eV[i].value for i in inst.obV]

data1 = {'V': V, 'PD': PD, 'QD': QD, 'PS': PS, 'QS': QS}
df1 = pd.DataFrame(data1).map(format_value)

data2 = {'I': I, 'P': P, 'Q': Q}
df2 = pd.DataFrame(data2, 
                   index=pd.MultiIndex.from_tuples(ind)).map(format_value)

print(df1)
print()
print(df2)
print()
print("eV [*] := ")
df3 = []
for i in inst.obV:
    row = [i, inst.eV[i].value]
    df3.append(row)
        
columns = ['i', 'eV']
df3 = pd.DataFrame(df3, columns=columns)
print(df3)

print()
print("eI [*] := ")
df4 = []
for i, j in inst.olI: 
    row = [i, j, inst.eI[i, j].value]
    df4.append(row)
        
columns = ['i', 'j', 'eI']
df4 = pd.DataFrame(df4, columns=columns)
print(df4)

print()
df5 = []
for i, j in inst.olP: 
    row = [i, j, inst.eP[i, j].value, inst.eQ[i, j].value]
    df5.append(row)
        
columns = ['i', 'j', 'eP', 'eQ']
df5 = pd.DataFrame(df5, columns=columns)
print(df5)

print()
df6 = []
for i in inst.obPS: 
    row = [i, inst.ePS[i].value, inst.eQS[i].value]
    df6.append(row)
        
columns = ['i', 'ePS', 'eQS']
df6 = pd.DataFrame(df6, columns=columns)
print(df6)

print()
df7 = []
for i in inst.obPD: 
    row = [i, inst.ePD[i].value, inst.eQD[i].value]
    df7.append(row)
        
columns = ['i', 'ePD', 'eQD']
df7 = pd.DataFrame(df7, columns=columns)
print(df7)

print()
print('Minimos Quadrados: %.4f' % inst.minimos_quadrados())