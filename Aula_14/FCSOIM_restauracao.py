from pyomo.environ import *
import pandas as pd
import numpy as np
import time

start = time.time()

m = AbstractModel()

# declaração de conjuntos e parâmetros
m.ob  = Set(dimen=1)                   # conjunto de nós
m.ol  = Set(dimen=2)                   # conjunto de circuitos      
m.och = Set(dimen=2)                   # conjunto das chaves
m.od  = Set(dimen=1)                   # conjunto de níveis de demanda
m.oz  = Set(dimen=1)                   # conjunto de zonas 

m.obout = Set(m.ob, within=m.ob, initialize=[])
m.obin  = Set(m.ob, within=m.ob, initialize=[])

m.obOUT = Set(m.ob, within=m.ob, initialize=[])
m.obIN  = Set(m.ob, within=m.ob, initialize=[])

def Populate_In_and_Out_Ol(m):
    for i, j in m.ol:
        m.obin[j].add(i)
        m.obout[i].add(j)
m.Populate_In_and_Out_Ol = BuildAction(rule=Populate_In_and_Out_Ol)

def Populate_In_and_Out_Och(m):
    for i, j in m.och:
        m.obIN[j].add(i)
        m.obOUT[i].add(j)
m.Populate_In_and_Out_Och = BuildAction(rule=Populate_In_and_Out_Och)


m.Tb     = Param(m.ob)                      # tipo de barra 0: carga, 1: SE
m.Smax   = Param(m.ob)                      # potencia aparente máxima
m.zb     = Param(m.ob)                      # zona onde pertence o nó i
m.PD     = Param(m.ob, m.od)                # potência ativa de demanda no nó i
m.QD     = Param(m.ob, m.od)                # potência reativa de demanda no nó i
m.R      = Param(m.ol, mutable=True)        # resistência no circuito ij
m.X      = Param(m.ol, mutable=True)        # reatância no circuito ij
m.Z2     = Param(m.ol, mutable=True)        # impedância no circuito ij ao quadrado
m.Imax   = Param(m.ol)                      # limite máximo da magnitude de corrente no circuito ij
m.zl     = Param(m.ol)                      # zona onde pertence o circuito ij
m.alpha  = Param(m.od)                      # número de horas no nível d    
m.cls    = Param(m.od)                      # custo das perdas de energia no nível d
m.cr     = Param(m.od)                      # custo de corte de carga
m.Rch    = Param(m.och, mutable=True)       # resistência na chave ij    
m.Ichmax = Param(m.och)                     # limite máximo da magnitude de corrente na chave ij
m.ini    = Param(m.och)                     # estado inicial da chave ij
m.zener  = Param(m.oz)                      # estado final de operação da zona
m.Vnom   = Param()                          # magnitude da tensão nominal
m.Vmin   = Param(initialize=0.95*m.Vnom)    # magnitude de tensão mínima
m.Vmax   = Param(initialize=1.00*m.Vnom)    # magnitude de tensão máxima
m.cch    = Param(initialize=0.2)            # custo de chaveamento


# declaração das variáveis
m.Vqdr   = Var(m.ob, m.od)                  # variável que representta o quadrado de V[i]
m.PS     = Var(m.ob, m.od)                  # potência ativa fornecida pela subestação no nó i
m.QS     = Var(m.ob, m.od)                  # potência reativa forneda pela subestação no nó i
m.Iqdr   = Var(m.ol, m.od)                  # variável que representa o quadrado de I[i,j]
m.P      = Var(m.ol, m.od)                  # fluxo de potência ativa no circuito ij
m.Q      = Var(m.ol, m.od)                  # fluxo de potência reativa no circuito ij
m.Pch    = Var(m.och, m.od)                 # fluxo de potência ativa na chave ij    
m.Qch    = Var(m.och, m.od)                 # fluxo de potência ativa na chave ij
m.Ichqdr = Var(m.och, m.od)                 # variável que representa o quadrado de Ich[i,j]
m.w      = Var(m.och, domain=Binary)        # estado de operação da chave
m.x      = Var(m.oz, domain=Binary)         # estado de operação da zona


# definição da função objetivo
# minimização das perdas de energia
def perdas(m):
    return sum(m.cr[d]*m.alpha[d]*(sum(m.PD[i, d]*(1 - m.x[m.zb[i]]) for i in m.ob)) for d in m.od) + \
           m.cch * sum(m.w[i, j] for i, j in m.och if m.ini[i, j] == 0) + \
           m.cch * sum(1 - m.w[i, j] for i, j in m.och if m.ini[i, j] == 1) + \
           sum(m.cls[d]*m.alpha[d]*(sum(m.R[i, j]*m.Iqdr[i, j, d] for i, j in m.ol) + \
                                    sum(m.Rch[i, j] * m.Ichqdr[i, j, d] for i, j in m.och)) for d in m.od)
m.perdas = Objective(rule=perdas, sense=minimize)

# definição das restrições
def balance_potencia_ativa(m, i, d):
    return sum(m.P[j, i, d] for j in m.obin[i]) - sum(m.P[i,j,d] + m.R[i,j] * m.Iqdr[i,j,d] for j in m.obout[i]) + \
           sum(m.Pch[j, i, d] for j in m.obIN[i]) - sum(m.Pch[i, j, d] for j in m.obOUT[i]) + m.PS[i, d] == m.PD[i, d] * m.x[m.zb[i]]
m.balance_potencia_ativa = Constraint(m.ob, m.od, rule=balance_potencia_ativa)

def balance_potencia_reativa(m, i, d):
    return sum(m.Q[j, i, d] for j in m.obin[i]) - sum(m.Q[i,j,d] + m.X[i,j] * m.Iqdr[i,j,d] for j in m.obout[i]) + \
           sum(m.Qch[j, i, d] for j in m.obIN[i]) - sum(m.Qch[i, j, d] for j in m.obOUT[i]) + m.QS[i, d] == m.QD[i, d] * m.x[m.zb[i]]
m.balance_potencia_reativa = Constraint(m.ob, m.od, rule=balance_potencia_reativa)

def queda_magnitude_tensao(m, i, j, d):
    return m.Vqdr[i,d] - 2*(m.R[i,j]*m.P[i,j,d] + m.X[i,j]*m.Q[i,j,d]) - m.Z2[i,j]*m.Iqdr[i,j,d] - m.Vqdr[j,d] == 0
m.queda_magnitude_tensao = Constraint(m.ol, m.od, rule=queda_magnitude_tensao)

def calculo_magnitude_corrente(m, i, j, d):
    return m.Vqdr[j,d] * m.Iqdr[i,j,d] >= m.P[i,j,d]**2 + m.Q[i,j,d]**2
m.calculo_magnitude_corrente = Constraint(m.ol, m.od, rule=calculo_magnitude_corrente)

def calculo_magnitude_corrente_chave(m, i, j, d):
    return m.Vqdr[j,d] * m.Ichqdr[i,j,d] >= m.Pch[i,j,d]**2 + m.Qch[i,j,d]**2
m.calculo_magnitude_corrente_chave = Constraint(m.och, m.od, rule=calculo_magnitude_corrente_chave)

def operacao_chave_1(m, i, j, d):
    return m.Vqdr[j,d] - m.Vqdr[i,d] <= 2 * m.Vmax**2 * (1 - m.w[i,j])
m.operacao_chave_1 = Constraint(m.och, m.od, rule=operacao_chave_1)

def operacao_chave_2(m, i, j, d):
    return - 2 * m.Vmax**2 * (1 - m.w[i,j]) <= m.Vqdr[j,d] - m.Vqdr[i,d]
m.operacao_chave_2 = Constraint(m.och, m.od, rule=operacao_chave_2)

def operacao_chave_3(m, i, j, d):
    return m.Pch[i,j,d] <= (m.Vmax * m.Ichmax[i,j]) * m.w[i,j]
m.operacao_chave_3 = Constraint(m.och, m.od, rule=operacao_chave_3)

def operacao_chave_4(m, i, j, d):
    return - (m.Vmax * m.Ichmax[i,j]) * m.w[i,j] <= m.Pch[i,j,d]
m.operacao_chave_4 = Constraint(m.och, m.od, rule=operacao_chave_4)

def operacao_chave_5(m, i, j, d):
    return m.Qch[i,j,d] <= (m.Vmax * m.Ichmax[i,j]) * m.w[i,j]
m.operacao_chave_5 = Constraint(m.och, m.od, rule=operacao_chave_5)

def operacao_chave_6(m, i, j, d):
    return - (m.Vmax * m.Ichmax[i,j]) * m.w[i,j] <= m.Qch[i,j,d]
m.operacao_chave_6 = Constraint(m.och, m.od, rule=operacao_chave_6)

def condicao_1_radialidade(m):
    return len(m.ol) + sum(m.w[i,j] for i, j in m.och) == len(m.ob) - sum(1 for i in m.ob if m.Tb[i]) - \
                                                                    - sum(1 - m.x[z] for z in m.oz)
m.condicao_1_radialidade = Constraint(rule=condicao_1_radialidade)

def limite_magnitude_tensão_1(m, i, d):
    return m.Vqdr[i,d] <= m.Vmax**2 * m.x[m.zb[i]]
m.limite_magnitude_tensão_1 = Constraint(m.ob, m.od, rule=limite_magnitude_tensão_1)

def limite_magnitude_tensão_2(m, i, d):
    return m.x[m.zb[i]] * m.Vmin**2 <= m.Vqdr[i,d]
m.limite_magnitude_tensão_2 = Constraint(m.ob, m.od, rule=limite_magnitude_tensão_2)

def limite_magnitude_corrente_1(m, i, j, d):
    return m.Iqdr[i,j,d] <= m.Imax[i,j]**2 * m.x[m.zl[i,j]]
m.limite_magnitude_corrente_1 = Constraint(m.ol, m.od, rule=limite_magnitude_corrente_1)

def limite_magnitude_corrente_2(m, i, j, d):
    return 0 <= m.Iqdr[i,j,d]
m.limite_magnitude_corrente_2 = Constraint(m.ol, m.od, rule=limite_magnitude_corrente_2)

def limite_magnitude_corrente_chaves_1(m, i, j, d):
    return m.Ichqdr[i,j,d] <= (m.Ichmax[i,j]**2) * m.w[i,j]
m.limite_magnitude_corrente_chaves_1 = Constraint(m.och, m.od, rule=limite_magnitude_corrente_chaves_1)

def limite_magnitude_corrente_chaves_2(m, i, j, d):
    return 0 <= m.Ichqdr[i,j,d]
m.limite_magnitude_corrente_chaves_2 = Constraint(m.och, m.od, rule=limite_magnitude_corrente_chaves_2)

def capacidade_subestacao(m, i, d):
    if m.Tb[i] == 1:
        return m.PS[i,d]**2 + m.QS[i,d]**2 <= m.Smax[i]**2
    else:
        return Constraint.Skip
m.capacidade_subestacao = Constraint(m.ob, m.od, rule=capacidade_subestacao)

def condicao_operacao_chave_zona_1(m, z):
    return sum(m.w[i,j] for i, j in m.och if m.zb[i] == z) + \
           sum(m.w[k,i] for k, i in m.och if m.zb[i] == z) >= m.x[z]
m.condicao_operacao_chave_zona_1 = Constraint(m.oz, rule=condicao_operacao_chave_zona_1)

def condicao_operacao_chave_zona_2(m, i, j):
    return m.w[i,j] <= m. x[m.zb[i]]
m.condicao_operacao_chave_zona_2 = Constraint(m.och, rule=condicao_operacao_chave_zona_2)

def condicao_operacao_chave_zona_3(m, i, j):
    return m.w[i,j] <= m.x[m.zb[j]]
m.condicao_operacao_chave_zona_3 = Constraint(m.och, rule=condicao_operacao_chave_zona_3)

def zona_energizada(m, z):
    return m.x[z] <= m.zener[z]
m.zona_energizada = Constraint(m.oz, rule=zona_energizada)


file = 'data/sistema44nos.dat'
inst = m.create_instance(file)


# fixar a magnitude de tensão da subestação no valor nominal
for i in inst.ob:
    for d in inst.od:
        if inst.Tb[i] == 1:
            inst.Vqdr[i, d].fix(inst.Vnom**2)
            
# fixar a injeção de potencia ativa nos nós de carga iguais a zero
for i in inst.ob:
    for d in inst.od:
        if inst.Tb[i] == 0:
            inst.PS[i, d].fix(0)
            
# fixar a injeção de potencia reativa nos nós de carga iguais a zero
for i in inst.ob:
    for d in inst.od:
        if inst.Tb[i] == 0:
            inst.QS[i, d].fix(0)        

# fixar a injeção de potencia reativa nos nós de carga iguais a zero
for z in inst.oz:
    if z == 0:
        inst.x[z].fix(1)
        
# nomalização e calculo da impedancia ao quadrado
for i, j in inst.ol:
    inst.R[i, j] = inst.R[i, j] / 1000      # de ohm para kohm
    inst.X[i, j] = inst.X[i, j] / 1000      # de ohm para kohm
    inst.Z2[i, j] = inst.R[i, j]**2 + inst.X[i, j]**2

for i, j in inst.och:
    inst.Rch[i, j] = inst.Rch[i, j] / 1000  # de ohm para kohm
    

solver = SolverFactory('cplex')
solver.options['mipgap'] = 0.0
results = solver.solve(inst, tee=True)
end = time.time()
print('Perdas: %.4f' % inst.perdas())


print("\nRESUMO\n")
print("Nivel de demanda                      :", end='')
for d in inst.od:
    print(f'{d: >10d}', end='')
print()
print("Perdas de potencia ativa (kW)         :", end='')
for d in inst.od:
    loss_p = sum(inst.R[i,j].value * inst.Iqdr[i,j,d].value for i,j in inst.ol)
    print(f'{loss_p: >10.2f}', end='')
print()
print("Perdas de potencia reativa (kVAr)     :", end='')
for d in inst.od:
    loss_q = sum(inst.X[i,j].value * inst.Iqdr[i,j,d].value for i,j in inst.ol)
    print(f'{loss_q: >10.2f}', end='')
print()
print("Magnitude de tensao minima (kV)       :", end='')
for d in inst.od:
    min_V = min(inst.Vqdr[j, d].value for j in inst.ob)
    print(f'{min_V ** 0.5: >10.4f}', end='')
print()
print("No da magnitude de tensao minima      :", end='')
for d in inst.od:
    min_i = min(i for i in inst.ob if inst.Vqdr[i, d].value == min(inst.Vqdr[j, d].value for j in inst.ob))
    print(f'{min_i: >10d}', end='')
print()
print("Potencia ativa da subestacao (kW)     :", end='')
for d in inst.od:
    min_i = min(i for i in inst.ob if inst.Tb[i] == 1)
    print(f'{inst.PS[min_i, d].value: >10.2f}', end='')
print()
print("Potencia reativa da subestacao (kVAr) :", end='')
for d in inst.od:
    min_i = min(i for i in inst.ob if inst.Tb[i] == 1)
    print(f'{inst.QS[min_i, d].value: >10.2f}', end='')
print()
print("Demanda de potencia ativa (kW)        :", end='')
for d in inst.od:
    demand = sum(inst.PD[i, d] for i in inst.ob)
    print(f'{demand: >10.2f}', end='')
print()
print("Demanda de potencia reativa (kVAr)    :", end='')
for d in inst.od:
    demand = sum(inst.QD[i, d] for i in inst.ob)
    print(f'{demand: >10.2f}', end='')
print()

print()
print()
print("MAGNITUDES DE TENSAO (kV)")
print("    Nivel de demanda\n  No", end='')
for d in inst.od:
    print(f' {d: >8d}', end='')
print()
for i in inst.ob:
    print(f'{i: >4d}', end='')
    for d in inst.od:
        voltage = inst.Vqdr[i, d].value ** 0.5
        print(f' {voltage: >8.4f}', end='')
    print()

print()
print()
print("MAGNITUDES DE FLUXO DE CORRENTE (A)")
print("         Nivel de demanda\n Circuito", end='')
for d in inst.od:
    print(f' {d: >8d}', end='')
print()
for i, j in inst.ol:
    print(f'{i: >4d} {j: >4d}', end='')
    for d in inst.od:
        current_flow = inst.Iqdr[i, j, d].value ** 0.5
        print(f' {current_flow: >8.2f}', end='')
    print()
    
print()
print()
print("FLUXO DE POTENCIA ATIVA (kW)")
print("         Nivel de demanda\n Circuito", end='')
for d in inst.od:
    print(f' {d: >8d}', end='')
print()
for i, j in inst.ol:
    print(f'{i: >4d} {j: >4d}', end='')
    for d in inst.od:
        active_power = inst.P[i, j, d].value
        print(f' {active_power: >8.2f}', end='')
    print()    
    
print()  
print()
print("FLUXO DE POTENCIA REATIVA (kVAr)")
print("         Nivel de demanda\n Circuito", end='')
for d in inst.od:
    print(f' {d: >8d}', end='')
print()
for i, j in inst.ol:
    print(f'{i: >4d} {j: >4d}', end='')
    for d in inst.od:
        reactive_power = inst.Q[i, j, d].value
        print(f' {reactive_power: >8.2f}', end='')
    print() 
    
print()
print("PERDAS DE POTENCIA ATIVA (kW)")
print("         Nivel de demanda\n Circuito", end='')
for d in inst.od:
    print(f' {d: >8d}', end='')
print()
for i, j in inst.ol:
    print(f'{i: >4d} {j: >4d}', end='')
    for d in inst.od:
        R =  inst.R[i, j].value
        Iqdr = inst.Iqdr[i, j, d].value
        active_power_loss = R * Iqdr 
        print(f' {active_power_loss: >8.2f}', end='')
    print()
    
print()
print(f"\nElapsed time:  {end-start:10.2f}s\n")
print()
df1 = []
for i, j in inst.och:
    row = [i, j, inst.ini[i,j], int(inst.w[i,j].value)]
    df1.append(row)
df1 = pd.DataFrame(df1, columns=['i', 'j', 'ini', 'w'])
print(df1)
    
print()
print('x[*]  :=  ')
for i in inst.oz:
    print(i, '  ', int(inst.x[i].value))
    