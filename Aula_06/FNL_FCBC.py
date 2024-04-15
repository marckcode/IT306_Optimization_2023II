from pyomo.environ import *
from time import time
t1 = time()

m = AbstractModel()

# declaração de conjuntos e parâmetros
m.ob = Set()                     # conjunto de nós
m.ol = Set(dimen=2)              # conjunto de circuitos
m.od = Set()                     # conjunto de níveis de demanda

m.obout = Set(m.ob, within=m.ob, initialize=[])
m.obin = Set(m.ob, within=m.ob, initialize=[])

def Populate_In_and_Out(m):
    for i, j in m.ol:
        m.obin[j].add(i)
        m.obout[i].add(j)
m.In_n_Out = BuildAction(rule=Populate_In_and_Out)


m.Tb    = Param(m.ob)                       # tipo de barra 0: carga, 1: SE
m.bsh   = Param(m.ob)                       # susceptância shunt
m.PD    = Param(m.ob, m.od)                 # potência ativa de demanda no nó i
m.QD    = Param(m.ob, m.od)                 # potência reativa de demanda no nó i
m.R     = Param(m.ol, mutable=True)         # resistência no circuito ij
m.X     = Param(m.ol, mutable=True)         # reatância no circuito ij
m.Imax  = Param(m.ol)                       # limite máximo da magnitude de corrente no circuito ij
m.alpha = Param(m.od)                       # número de horas no nível d
m.cls   = Param(m.od)                       # custo das perdas de energia no nível d
m.Vnom  = Param(mutable=True, initialize=11.0/sqrt(3))      # magnitude da tensão nominal

# declaração das variáveis
m.Vre  = Var(m.ob, m.od)            # variável que representa o quadrado de V[i,d]
m.Vim  = Var(m.ob, m.od)            # potência ativa fornecida pela subestação no nó i
m.ISre = Var(m.ob, m.od)            # variável que representa o quadrado de I[i,j,d]  
m.ISim = Var(m.ob, m.od)            # potência reativa fornecida pela subestação no nó i
m.IDre = Var(m.ob, m.od)            # fluxo de potência ativa no circuito ij
m.IDim = Var(m.ob, m.od)            # fluxo de potência reativa no circuito ij
m.Ire  = Var(m.ol, m.od)            # Parte real do fluxo da corrente 
m.Iim  = Var(m.ol, m.od)            # Parte imaginaria do fluxo da corrente    

# definição da função objetivo
def custo(m):
    return (sum(m.cls[d]*m.alpha[d]*
           sum(m.Vre[i,d]*m.ISre[i,d] + m.Vim[i,d]*m.ISim[i,d] 
                        for i in m.ob if m.Tb[i]==1) for d in m.od)) 
m.custo = Objective(rule=custo, sense=minimize)

def balance_corrente_real(m, i, d):
    return sum(m.Ire[j,i,d] for j in m.obin[i]) - sum(m.Ire[i,j,d] for j in m.obout[i]) + m.ISre[i,d] == m.IDre[i,d]
m.balance_corrente_real = Constraint(m.ob, m.od, rule=balance_corrente_real)

def balance_corrente_imaginaria(m, i, d):
    return sum(m.Iim[j,i,d] for j in m.obin[i]) - sum(m.Iim[i,j,d] for j in m.obout[i]) + m.ISim[i,d] == m.IDim[i,d]
m.balance_corrente_imaginaria = Constraint(m.ob, m.od, rule=balance_corrente_imaginaria)

def queda_tensao_real(m, i, j, d):
    return m.Vre[i,d] - m.Vre[j,d] == m.Ire[i,j,d] * m.R[i,j] - m.Iim[i,j,d] * m.X[i,j]
m.queda_tensao_real = Constraint(m.ol, m.od, rule=queda_tensao_real)

def queda_tensao_imaginaria(m, i, j, d):
    return m.Vim[i,d] - m.Vim[j,d] == m.Ire[i,j,d] * m.X[i,j] + m.Iim[i,j,d] * m.R[i,j]
m.queda_tensao_imaginaria = Constraint(m.ol, m.od, rule=queda_tensao_imaginaria)

def demanda_potencia_ativa(m, i, d):
    return m.PD[i, d] == m.Vre[i, d] * m.IDre[i, d] + m.Vim[i, d] * m.IDim[i, d]
m.demanda_potencia_ativa = Constraint(m.ob, m.od, rule=demanda_potencia_ativa)

def demanda_potencia_reativa(m, i, d):
    return m.QD[i, d] == -m.Vre[i, d] * m.IDim[i, d] + m.Vim[i, d] * m.IDre[i, d]
m.demanda_potencia_reativa = Constraint(m.ob, m.od, rule=demanda_potencia_reativa)


file = "data/sistema34nos.dat"
instance = m.create_instance(file)

# fixar a parte real da tensão da subestação no valor nominal
for i in instance.ob:
    for d in instance.od:
        if instance.Tb[i] == 1:
            instance.Vre[i, d].fix(instance.Vnom)
            
# fixar a parte imaginaria da tensão da subestação a zero
for i in instance.ob:
    for d in instance.od:
        if instance.Tb[i] == 1:
            instance.Vim[i, d].fix(0)
            
# fixar a injeção de potencia ativa nos nós de carga iguais a zero
for i in instance.ob:
    for d in instance.od:
        if instance.Tb[i] == 0:
            instance.ISre[i, d].fix(0)     
            
# fixar a injeção de potencia reativa nos nós de carga iguais a zero
for i in instance.ob:
    for d in instance.od:
        if instance.Tb[i] == 0:
            instance.ISim[i, d].fix(0)    
            
# nomalização e calculo da impedancia ao quadrado
for i, j in instance.ol:
    instance.R[i, j] = instance.R[i, j] / 1000      # de ohm para kohm
    instance.X[i, j] = instance.X[i, j] / 1000      # de ohm para kohm


solver = SolverFactory('ipopt')
results = solver.solve(instance, tee=True)
print('Custo total: %.4f' % instance.custo())

print("\nRESUMO")
print("Nivel de demanda                      :", end='')
for d in instance.od:
    print(f'{d: >10d}', end='')
print()

print("Perdas de potencia ativa (kW)         :", end='')
for d in instance.od:
    loss = sum(instance.R[i,j].value * (instance.Ire[i,j,d].value**2 + instance.Iim[i,j,d].value**2) for (i,j) in instance.ol)
    print(f'{loss: >10.2f}', end='')
print()

print("Perdas de potencia reativa (kVAr)     :", end='')
for d in instance.od:
    loss = sum(instance.X[i,j].value * (instance.Ire[i,j,d].value**2 + instance.Iim[i,j,d].value**2) for (i,j) in instance.ol)
    print(f'{loss: >10.2f}', end='')
print()

print("Magnitude de tensao minima (kV)       :", end='')
for d in instance.od:
    min_V = min(instance.Vre[j,d].value**2+instance.Vim[j,d].value**2 for j in instance.ob)
    print(f'{min_V ** 0.5: >10.4f}', end='')
print()

print("Magnitude de tensao minima (pu)       :", end='')
for d in instance.od:
    min_V = min(instance.Vre[j,d].value**2+instance.Vim[j,d].value**2 for j in instance.ob)
    print(f'{(min_V ** 0.5) / instance.Vnom.value: >10.4f}', end='')
print()

print("No da magnitude de tensao minima      :", end='')
for d in instance.od:
    min_i = min(i for i in instance.ob if instance.Vre[i,d].value**2+instance.Vim[i,d].value**2 == 
            min(instance.Vre[j, d].value**2+instance.Vim[j,d].value**2 for j in instance.ob))
    print(f'{min_i: >10d}', end='')
print()

print("Potencia ativa da subestacao (kW)     :", end='')
for d in instance.od:
    ind = sum(instance.Vre[i,d].value*instance.ISre[i,d].value 
              + instance.Vim[i,d].value*instance.ISim[i,d].value for i in instance.ob if instance.Tb[i] == 1)
    print(f'{ind: >10.2f}', end='')
print()

print("Potencia reativa da subestacao (kVAr) :", end='')
for d in instance.od:
    ind = sum(-instance.Vre[i,d].value*instance.ISim[i,d].value 
              + instance.Vim[i,d].value*instance.ISre[i,d].value for i in instance.ob if instance.Tb[i] == 1)
    print(f'{ind: >10.2f}', end='')
print()

print("Demanda de potencia ativa (kW)        :", end='')
for d in instance.od:
    demand = sum(instance.PD[i, d] for i in instance.ob)
    print(f'{demand: >10.2f}', end='')
print()

print("Demanda de potencia reativa (kVAr)    :", end='')
for d in instance.od:
    demand = sum(instance.QD[i, d] for i in instance.ob)
    print(f'{demand: >10.2f}', end='')
print()


print()
print()
print("MAGNITUDES DE TENSAO (kV)")
print("    Nivel de demanda\n  No", end='')
for d in instance.od:
    print(f' {d: >8d}', end='')
print()
for i in instance.ob:
    print(f'{i: >4d}', end='')
    for d in instance.od:
        voltage = (instance.Vre[i,d].value**2 + instance.Vim[i,d].value**2)** 0.5
        print(f' {voltage: >8.4f}', end='')
    print()
    
    
print()
print()
print("MAGNITUDES DE FLUXO DE CORRENTE (A)")
print("         Nivel de demanda\n Circuito", end='')
for d in instance.od:
    print(f' {d: >8d}', end='')
print("     Imax\n")
for i, j in instance.ol:
    print(f'{i: >4d} {j: >4d}', end='')
    for d in instance.od:
        current_flow = (instance.Ire[i,j,d].value**2 + instance.Iim[i,j,d].value**2) ** 0.5
        print(f' {current_flow: >8.2f}', end='')
    imax_value = instance.Imax[i, j]
    print(f' {imax_value: >8.2f}')
    
    
print()
print()
print("FLUXO DE POTENCIA ATIVA (kW)")
print("         Nivel de demanda\n Circuito", end='')
for d in instance.od:
    print(f' {d: >8d}', end='')
print()
for i, j in instance.ol:
    print(f'{i: >4d} {j: >4d}', end='')
    for d in instance.od:
        active_power = instance.Vre[j,d].value * instance.Ire[i,j,d].value \
                       + instance.Vim[j,d].value*instance.Iim[i,j,d].value
        print(f' {active_power: >8.2f}', end='')
    print() 
    
    
print()  
print()
print("FLUXO DE POTENCIA REATIVA (kVAr)")
print("         Nivel de demanda\n Circuito", end='')
for d in instance.od:
    print(f' {d: >8d}', end='')
print()

# Print reactive power for each (i, j) and d
for i, j in instance.ol:
    print(f'{i: >4d} {j: >4d}', end='')
    for d in instance.od:
        reactive_power = -instance.Vre[j,d].value * instance.Iim[i,j,d].value \
                         + instance.Vim[j,d].value*instance.Ire[i,j,d].value
        print(f' {reactive_power: >8.2f}', end='')
    print()    
    

print()
print("PERDAS DE POTENCIA ATIVA (kW)")
print("         Nivel de demanda\n Circuito", end='')
for d in instance.od:
    print(f' {d: >8d}', end='')
print()

# Print active power losses for each (i, j) and d
for i, j in instance.ol:
    print(f'{i: >4d} {j: >4d}', end='')
    for d in instance.od:
        active_power_loss = instance.R[i,j].value * (instance.Ire[i,j,d].value**2 + instance.Iim[i,j,d].value**2) 
        print(f' {active_power_loss: >8.2f}', end='')
    print()

t2 = time()

print()
print()
print(f"Elapsed time  : {(t2 - t1): >10.2f}s")
print()
print()








