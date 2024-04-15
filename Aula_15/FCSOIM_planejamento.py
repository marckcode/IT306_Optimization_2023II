from pyomo.environ import *
import pandas as pd
import time

start = time.time()

m = AbstractModel()

m.Ob = Set()                    # set of nodes
m.Ol = Set(dimen=2)             # set of branches 
m.Oa = Set()                    # set of types of conductors

m.obout = Set(m.Ob, within=m.Ob, initialize=[])
m.obin  = Set(m.Ob, within=m.Ob, initialize=[])

def Populate_In_and_Out_Ol(m):
    for i, j in m.Ol:
        m.obin[j].add(i)
        m.obout[i].add(j)
m.Populate_In_and_Out_Ol = BuildAction(rule=Populate_In_and_Out_Ol)

m.Tp   = Param(m.Ob)                        # node type
m.Sd   = Param(m.Ob)                        # apparent demand
m.Pd   = Param(m.Ob, mutable=True)          # active demand
m.Qd   = Param(m.Ob, mutable=True)          # reactive demand 
m.S0   = Param(m.Ob)                        # existent apparent power 
m.SL   = Param(m.Ob)                        # capacity of construction or increasing of apparent power
m.CCS  = Param(m.Ob)                        # construction or increasing of apparent power
m.RT   = Param(m.Oa, mutable=True)          # resistance (ohm/km)
m.XT   = Param(m.Oa, mutable=True)          # reactance (ohm/km)
m.ZT2  = Param(m.Oa, mutable=True)          # impedance (ohm/km)
m.IT   = Param(m.Oa)                        # max current (A)
m.CFT  = Param(m.Oa)                        # cost of conductor ($/km)
m.Di   = Param(m.Ol)                        # length of the branch (km)
m.Vnom = Param(initialize=34.5)             # nominal voltage magnitude (kV)
m.Vmax = Param(initialize=1.03 * m.Vnom)    # upper bound of voltage (kV)
m.Vmin = Param(initialize=0.97 * m.Vnom)    # lower bound of voltage (kV)
m.pf   = Param(initialize=0.9)              # power factor
m.LF   = Param(initialize=0.35)             # load factor
m.cl   = Param(initialize=0.05)             # cost of energy (US$/kWh)
m.cs   = Param(initialize=0.0)              # operation cost of substation (US$/kWh)
m.alfa = Param(initialize=8760)             # number of hours in one year
m.K    = Param(initialize=20)               # planning years;
m.tj   = Param(initialize=0.1)              # interest rate

m.V2  = Var(m.Ob, domain=NonNegativeReals)
m.Sg2 = Var(m.Ob, domain=NonNegativeReals)
m.Pg  = Var(m.Ob, domain=NonNegativeReals)
m.Qg  = Var(m.Ob, domain=NonNegativeReals)
m.I2s = Var(m.Ol, domain=NonNegativeReals)
m.Ps  = Var(m.Ol)
m.Qs  = Var(m.Ol)
m.b   = Var(m.Ol)
m.I2  = Var(m.Ol, m.Oa, domain=NonNegativeReals)
m.P   = Var(m.Ol, m.Oa)
m.Q   = Var(m.Ol, m.Oa)
m.z   = Var(m.Ol, m.Oa, domain=Binary)
m.xg  = Var(m.Ob, domain=Binary, initialize=0)

def investimento(m):
    return (sum(m.z[i, j, a] * m.CFT[a] * m.Di[i, j] for i, j in m.Ol for a in m.Oa) + 
            sum(m.CCS[i] * m.xg[i] for i in m.Ob if m.Tp[i] == 1) + (sum(m.alfa * m.LF * m.cl * m.RT[a] * 
                m.I2[i, j, a] * m.Di[i, j] for i, j in m.Ol for a in m.Oa) * (1 / m.tj) * (1 - (1 / ((1 + m.tj) ** m.K)))) +
           (sum(m.alfa * m.LF * m.cs * m.Sg2[i] for i in m.Ob if m.Tp[i] == 1) * (1 / m.tj) * (1 - (1 / ((1 + m.tj) ** m.K)))))
m.investimento = Objective(rule=investimento, sense=minimize)

def active_power_balance(m, i):
    return (- sum(m.P[i,j,a] + m.I2[i,j,a] * m.RT[a] * m.Di[i,j] for j in m.obout[i] for a in m.Oa) +
              sum(m.P[j,i,a] for j in m.obin[i] for a in m.Oa) + m.Pg[i] == m.Pd[i])
m.active_power_balance = Constraint(m.Ob, rule=active_power_balance)

def reactive_power_balance(m, i):
    return (- sum(m.Q[i,j,a] + m.I2[i,j,a] * m.XT[a] * m.Di[i,j] for j in m.obout[i] for a in m.Oa) +
              sum(m.Q[j,i,a] for j in m.obin[i] for a in m.Oa) + m.Qg[i] == m.Qd[i])
m.reactive_power_balance = Constraint(m.Ob, rule=reactive_power_balance)

def voltage_drop(m, i, j):
    return (m.V2[i] - m.V2[j] == m.b[i,j] + 
            sum(2 * (m.P[i,j,a] * m.RT[a] + m.Q[i,j,a] * m.XT[a]) * m.Di[i,j] + 
                m.I2[i,j,a] * m.ZT2[a] * m.Di[i,j] ** 2 for a in m.Oa))
m.voltage_drop = Constraint(m.Ol, rule=voltage_drop)

def current_flow_square(m, i, j):
    return m.V2[j] * m.I2s[i,j] >= m.Ps[i,j]**2 + m.Qs[i,j]**2
m.current_flow_square = Constraint(m.Ol, rule=current_flow_square)

def I_soma(m, i, j):
    return m.I2s[i,j] == sum(m.I2[i,j,a] for a in m.Oa)
m.I_soma = Constraint(m.Ol, rule=I_soma)

def P_soma(m, i, j):
    return m.Ps[i,j] == sum(m.P[i,j,a] for a in m.Oa)
m.P_soma = Constraint(m.Ol, rule=P_soma)

def Q_soma(m, i, j):
    return m.Qs[i,j] == sum(m.Q[i,j,a] for a in m.Oa)
m.Q_soma = Constraint(m.Ol, rule=Q_soma)

def VM_limits(m, i):
    return inequality(m.Vmin ** 2, m.V2[i], m.Vmax ** 2)
m.VM_limits = Constraint(m.Ob, rule=VM_limits)

def IM_max(m, i, j, a):
    return m.I2[i,j,a] <= m.IT[a] ** 2 * m.z[i,j,a]
m.IM_max = Constraint(m.Ol, m.Oa, rule=IM_max)

def limit_b_neg(m, i, j):
    return m.b[i,j] >= -(m.Vmax**2 - m.Vmin**2) * (1 - sum(m.z[i,j,a] for a in m.Oa))
m.limit_b_neg = Constraint(m.Ol, rule=limit_b_neg)

def limit_b_pos(m, i, j):
    return m.b[i,j] <= (m.Vmax**2 - m.Vmin**2) * (1 - sum(m.z[i,j,a] for a in m.Oa))
m.limit_b_pos = Constraint(m.Ol, rule=limit_b_pos)

def active_power_limit_z1(m, i, j, a):
    return m.P[i,j,a] <= m.Vmax * m.IT[a] * m.z[i,j,a]
m.active_power_limit_z1 = Constraint(m.Ol, m.Oa, rule=active_power_limit_z1)

def active_power_limit_z2(m, i, j, a):
    return m.P[i,j,a] >= - m.Vmax * m.IT[a] * m.z[i,j,a]
m.active_power_limit_z2 = Constraint(m.Ol, m.Oa, rule=active_power_limit_z2)

def reactive_power_limit_z1(m, i, j, a):
    return m.Q[i,j,a] <= m.Vmax * m.IT[a] * m.z[i,j,a]
m.reactive_power_limit_z1 = Constraint(m.Ol, m.Oa, rule=reactive_power_limit_z1)

def reactive_power_limit_z2(m, i, j, a):
    return m.Q[i,j,a] >= - m.Vmax * m.IT[a] * m.z[i,j,a]
m.reactive_power_limit_z2 = Constraint(m.Ol, m.Oa, rule=reactive_power_limit_z2)

def SE_capacity1(m, i):
    if m.Tp[i] == 1:
        return m.Sg2[i] >= m.Pg[i]**2 + m.Qg[i]**2
    else:
        return Constraint.Skip
m.SE_capacity1 = Constraint(m.Ob, rule=SE_capacity1)

def SE_capacity_total(m, i):
    if m.Tp[i] == 1:
        return m.Sg2[i] <= (m.S0[i] + m.SL[i] * m.xg[i]) ** 2
    else:
        return Constraint.Skip
m.SE_capacity_total = Constraint(m.Ob, rule=SE_capacity_total)

def uniqueness(m, i, j):
    return sum(m.z[i,j,a] for a in m.Oa) <= 1
m.uniqueness = Constraint(m.Ol, rule=uniqueness)

def three_network(m):
    return sum(m.z[i,j,a] for i, j in m.Ol for a in m.Oa) == len(m.Ob) - sum(1 for i in m.Ob if m.Tp[i] == 1)
m.three_network = Constraint(rule=three_network)


file = 'data/sistema23nos.dat'
inst = m.create_instance(file)

for a in inst.Oa:
    inst.RT[a] = inst.RT[a] / 1000 # para kohm
    inst.XT[a] = inst.XT[a] / 1000 # para kohm
    inst.ZT2[a] = inst.RT[a]**2 + inst.XT[a]**2

for i in inst.Ob:
    inst.Pd[i] = inst.Sd[i] * inst.pf
    inst.Qd[i] = inst.Sd[i] * sin(acos(inst.pf))
    
    if inst.Tp[i] == 0:
        inst.Pg[i].fix(0)
        inst.Qg[i].fix(0)
        
        
solver = SolverFactory('cplex')
results = solver.solve(inst, tee=True)
end = time.time()
FO = inst.investimento()
print('Investimento: %.4f' % FO)   

# Impressao do resultados
CF_se  = sum(inst.z[i, j, a].value * inst.CFT[a] * inst.Di[i, j] for i, j in inst.Ol for a in inst.Oa)
CF_lin = sum(inst.CCS[i] * inst.xg[i].value for i in inst.Ob if inst.Tp[i] == 1)
C_CSE  = (sum(inst.alfa * inst.LF * inst.cs * inst.Sg2[i].value for i in inst.Ob if inst.Tp[i] == 1) * 
             (1 / inst.tj) * (1 - (1 / ((1 + inst.tj) ** inst.K))))
C_LOSS = (sum(inst.alfa * inst.LF * inst.cl * inst.RT[a].value * inst.I2[i, j, a].value * 
              inst.Di[i, j] for i, j in inst.Ol for a in inst.Oa) * (1 / inst.tj) * (1 - (1 / ((1 + inst.tj) ** inst.K))))
perda_ativa   = sum(inst.I2[m, n, a].value * inst.RT[a].value * inst.Di[m, n] for (m, n) in inst.Ol for a in inst.Oa)
perda_reativa = sum(inst.I2[m, n, a].value * inst.XT[a].value * inst.Di[m, n] for (m, n) in inst.Ol for a in inst.Oa)

with open("Result.out", "w") as rf:
    rf.write("RESUMO\n")
    rf.write(f"Custo de investimento total [U$]:   {FO:12.2f}\n")
    rf.write(f"Custo investimento subestacoes [U$]:{CF_lin:12.2f}\n")
    rf.write(f"Custo investimento circuitos [U$]:  {CF_se:12.2f}\n")
    rf.write(f"Custo operacao subestacoes [U$]:    {C_CSE:12.2f}\n")
    rf.write(f"Custo perda potencia ativa [U$]:    {C_LOSS:12.2f}\n")
    rf.write(f"Geracao de potencia aparente [kVA]: {sum(sqrt(inst.Sg2[n].value) for n in inst.Ob if inst.Tp[n] == 1):12.2f}\n")
    rf.write(f"Geracao de potencia ativa [kW]:     {sum(inst.Pg[n].value for n in inst.Ob if inst.Tp[n] == 1):12.2f}\n")
    rf.write(f"Geracao de potencia reativa [kVAr]: {sum(inst.Qg[n].value for n in inst.Ob if inst.Tp[n] == 1):12.2f}\n")
    rf.write(f"Demanda de potencia ativa [kW]:     {sum(inst.Pd[n].value for n in inst.Ob):12.2f}\n")
    rf.write(f"Demanda de potencia reativa [kVAr]: {sum(inst.Qd[n].value for n in inst.Ob):12.2f}\n")
    rf.write(f"Perda potencia ativa [kW]:          {perda_ativa:12.2f}\n")
    rf.write(f"Perda potencia reativa [kVAr]:      {perda_reativa:12.2f}\n")
    rf.write(f"Magnitude de tensao minima [pu]:    {min(sqrt(inst.V2[n].value) / inst.Vnom for n in inst.Ob):12.4f}\n")
    
    rf.write("\nMAGNITUDES DE TENSAO\n")
    rf.write("  No     [pu]\n")
    for n in inst.Ob:
        rf.write(f"{n:4d} {sqrt(inst.V2[n].value) / inst.Vnom:8.4f}\n")

    rf.write("\nOPERACAO DAS SUBESTACOES\n")
    rf.write("  No       Pg       Qg       Sg     Smax   xg\n")
    rf.write("         [kW]   [kVAr]    [kVA]    [kVA]\n")
    for n in inst.Ob:
        if inst.Tp[n] == 1:
            rf.write(f"{n:4d} {inst.Pg[n].value:8.2f} {inst.Qg[n].value:8.2f} {sqrt(inst.Sg2[n].value):8.2f} {(inst.S0[n] + inst.SL[n] * inst.xg[n].value):8.2f} {inst.xg[n].value:4d}\n")

    rf.write("\nRESULTADOS DOS RAMOS\n")
    rf.write("   m   n       Vm       Vn      Pmn      Qmn       Imn Cond\n")
    rf.write("             [pu]     [pu]     [kW]    [kVAr]      [A]\n")
    for (m, n) in inst.Ol:
        rf.write(f"{m:4d} {n:4d} {sqrt(inst.V2[m].value) / inst.Vnom:8.4f} {sqrt(inst.V2[n].value) / inst.Vnom:8.4f} ")
        rf.write(f"{sum(inst.P[m,n,a].value for a in inst.Oa):8.2f} {sum(inst.Q[m,n,a].value for a in inst.Oa):8.2f} ")
        rf.write(f"{sum(sqrt(inst.I2[m,n,a].value) for a in inst.Oa):8.2f} {int(sum(a * inst.z[m,n,a].value for a in inst.Oa)):4d}\n")

    rf.write(f"\nElapsed time:  {end-start:10.2f}s\n")
