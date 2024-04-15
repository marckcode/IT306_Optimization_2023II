from pyomo.environ import *
import pandas as pd
import numpy as np

m = AbstractModel()

# declaração de conjuntos e parâmetros
m.ob   = Set()                      # conjunto de nós
m.ol   = Set(dimen=2)               # conjunto de circuitos
m.obV  = Set()                      # conjunto de medidas de magnitude de tensão
m.obPS = Set()                      # conjunto de medidas de potencia ativa da subestação
m.obPD = Set()                      # conjunto de medidas de potencia ativa da demanda
m.obQD = Set()                      # conjunto de medidas de potencia reativa da demanda

m.obout = Set(m.ob, within=m.ob, initialize=[])
m.obin = Set(m.ob, within=m.ob, initialize=[])

def Populate_In_and_Out(m):
    for i, j in m.ol:
        m.obin[j].add(i)
        m.obout[i].add(j)
m.In_n_Out = BuildAction(rule=Populate_In_and_Out)

m.Tb    = Param(m.ob)                   # tipo de barra 0: carga, 1: SE
m.Tp    = Param(m.ob)                   # barra de passagem 1, caso contrario 0
m.R     = Param(m.ol, mutable=True)     # resistência no circuito ij
m.X     = Param(m.ol, mutable=True)     # reatância no circuito ij
m.Z2    = Param(m.ol, mutable=True)     # impedância no circuito ij ao quadrado
m.Imax  = Param(m.ol)                   # limite máximo da magnitude de corrente no circuito ij
m.Vnom  = Param(mutable=True)           # magnitude da tensão nominal    
m.Vmax  = Param(mutable=True)           # magnitude máxima da tensão
m.Vmin  = Param(mutable=True)           # magnitude mínima da tensão
m.Vmed  = Param(m.obV)                  # medidas de magnitude de tensão
m.PSmed = Param(m.obPS)                 # medidas de potencia ativa da subestação
m.PDmed = Param(m.obPD)                 # medidas de potencia ativa da demanda
m.QDmed = Param(m.obQD)                 # medidas de potencia reativa da demanda
m.WV    = Param(m.obV)                  # ponderações das medidas de magnitude de tensão
m.WPS   = Param(m.obPS)                 # ponderações das medidas de potencia ativa da subestação
m.WPD   = Param(m.obPD)                 # ponderações das medidas de potencia ativa da demanda
m.WQD   = Param(m.obQD)                 # ponderações das medidas de potencia reativa da demanda

# declaração das variáveis
m.Vqdr  = Var(m.ob)                         # variável que representa a magnitude de tensão
m.PD    = Var(m.ob)                         # potência ativa de demanda no nó i
m.QD    = Var(m.ob)                         # potência reativa de demanda no nó i
m.PS    = Var(m.ob)                         # potência ativa fornecida pela subestação no nó i
m.QS    = Var(m.ob)                         # potência reativa fornecida pela subestação no nó i
m.Iqdr  = Var(m.ol)                         # variável que representa o fluxo de corrente
m.P     = Var(m.ol)                         # fluxo de potência ativa no circuito ij
m.Q     = Var(m.ol)                         # fluxo de potência reativa no circuito ij
m.eV    = Var(m.obV,  initialize=0)         # resíduos da magnitude de tensão 
m.ePS   = Var(m.obPS, initialize=0)         # resíduos da potencia ativa da subestação
m.ePD   = Var(m.obPD, initialize=0)         # resíduos da demanda de potencia ativa
m.eQD   = Var(m.obQD, initialize=0)         # resíduos da demanda de potencia reativa
m.w     = Var(m.obPD, initialize=0, 
                domain=NonNegativeReals)       # porcentagem de demanda de potência ativa não contabilizada 


# definição da função objetivo
def minimos_quadrados(m):
    return (sum(m.eV[k]*m.WV[k]*m.eV[k] for k in m.obV) + sum(m.ePS[k]*m.WPS[k]*m.ePS[k] for k in m.obPS) + 
            sum(m.ePD[k]*m.WPD[k]*m.ePD[k] for k in m.obPD) + sum(m.eQD[k]*m.WQD[k]*m.eQD[k] for k in m.obQD))
m.minimos_quadrados = Objective(rule=minimos_quadrados, sense=minimize)

# Definição das restrições
def calculo_residuos_tensao(m, k):
    return m.Vmed[k]**2 == m.Vqdr[k] + 2*m.Vmed[k]*m.eV[k]
m.calculo_residuos_tensao = Constraint(m.obV, rule=calculo_residuos_tensao)

def calculo_residuos_injecao_ativa(m, k):
    return m.PSmed[k] == m.PS[k] + m.ePS[k]
m.calculo_residuos_injecao_ativa = Constraint(m.obPS, rule=calculo_residuos_injecao_ativa)

def calculo_residuos_demanda_ativa(m, k):
    return m.PDmed[k] * (1 + m.w[k]) == m.PD[k] + m.ePD[k]
m.calculo_residuos_demanda_ativa = Constraint(m.obPD, rule=calculo_residuos_demanda_ativa)

def calculo_residuos_demanda_reativa(m, k):
    return m.QDmed[k] == m.QD[k] + m.eQD[k]
m.calculo_residuos_demanda_reativa = Constraint(m.obQD, rule=calculo_residuos_demanda_reativa)

def balance_potencia_ativa(m, i):
    return sum(m.P[j, i] for j in m.obin[i]) \
           - sum(m.P[i, j] + m.R[i, j] * m.Iqdr[i, j] for j in m.obout[i]) + m.PS[i] == m.PD[i]
m.balance_potencia_ativa = Constraint(m.ob, rule=balance_potencia_ativa)

def balance_potencia_reativa(m, i):
    return sum(m.Q[j, i] for j in m.obin[i]) \
           - sum(m.Q[i, j] + m.R[i, j] * m.Iqdr[i, j] for j in m.obout[i]) + m.QS[i] == m.QD[i]
m.balance_potencia_reativa = Constraint(m.ob, rule=balance_potencia_reativa)

def queda_magnitude_tensao(m, i, j):
    return m.Vqdr[i] - 2*(m.R[i, j]*m.P[i, j] + m.X[i, j]*m.Q[i, j]) - m.Z2[i, j]*m.Iqdr[i, j] - m.Vqdr[j] == 0
m.queda_magnitude_tensao = Constraint(m.ol, rule=queda_magnitude_tensao)

def calculo_magnitude_corrente(m, i, j):
    return m.Vqdr[j] * m.Iqdr[i, j] >= m.P[i, j]**2 + m.Q[i, j]**2
m.calculo_magnitude_corrente = Constraint(m.ol, rule=calculo_magnitude_corrente)

def limite_magnitude_corrente(m, i, j):
    return inequality(0, m.Iqdr[i, j], m.Imax[i, j]**2)
m.limite_magnitude_corrente = Constraint(m.ol, rule=limite_magnitude_corrente)

def limite_magnitude_tensão(m, i):
    return inequality(m.Vmin**2, m.Vqdr[i], m.Vmax**2)
m.limite_magnitude_tensão = Constraint(m.ob, rule=limite_magnitude_tensão)

file = 'data/FCSO_EE.dat'
inst = m.create_instance(file)

# tensão de fase
inst.Vnom = 0.220 / sqrt(3) # kV
inst.Vmax = 0.220 / sqrt(3) # kV
inst.Vmin = 0.200 / sqrt(3) # kV
    
# normalização e calculo da impedância ao quadrado
for i, j in inst.ol:
    inst.R[i, j] = inst.R[i, j] / 1000      # de ohm para kohm     
    inst.X[i, j] = inst.X[i, j] / 1000      # de ohm para kohm  
    inst.Z2[i, j] = inst.R[i, j]**2 + inst.X[i, j]**2
     

solver = SolverFactory('gurobi') # tem alguns problemas usando "cplex"
results = solver.solve(inst, tee=True)
print('Minimos cuadrados: %.2f' % inst.minimos_quadrados())

# definir uma funcão de customização
def format_value(value):
    if abs(value) == 0 or abs(value) == 1:
        return "."
    else:
        return abs(value)

matriz_ob = []
for t in inst.ob:
    row = [inst.Vqdr[t].value, inst.PD[t].value, inst.QD[t].value,
           inst.PS[t].value, inst.QS[t].value]
    matriz_ob.append(row)

columns = ['Vqdr', 'PD', 'QD', 'PS', 'QS']
matriz_ob = pd.DataFrame(matriz_ob, columns=columns)
print(matriz_ob)

ev  = {t:inst.eV[t].value for t in inst.obV}
eps = {t:inst.ePS[t].value for t in inst.obPS}
epd = {t:inst.ePD[t].value for t in inst.obPD}
w   = {t:inst.w[t].value for t in inst.obPD}
eqd = {t:inst.eQD[t].value for t in inst.obQD}
print("\n")

matriz_vsd = pd.DataFrame([ev, eps, epd, w, eqd]).T
matriz_vsd.fillna(0, inplace=True)
matriz_vsd.columns = ['eV', 'ePS', 'ePD', 'w', 'wQD']
matriz_vsd = matriz_vsd.map(format_value)
print(matriz_vsd)


