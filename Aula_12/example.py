from pyomo.environ import *

m = ConcreteModel('example')

m.f12 = Var()
m.f13 = Var()
m.f14 = Var()
m.f23 = Var()
m.f25 = Var()
m.f34 = Var()
m.f35 = Var()
m.f45 = Var()

m.x12 = Var(domain=Binary)
m.x13 = Var(domain=Binary)
m.x14 = Var(domain=Binary)
m.x23 = Var(domain=Binary)
m.x25 = Var(domain=Binary)
m.x34 = Var(domain=Binary)
m.x35 = Var(domain=Binary)
m.x45 = Var(domain=Binary)

m.g = Var()

def objective_function(model):
    return (m.f12*m.f12 + m.f13*m.f13 + m.f14*m.f14 + m.f23*m.f23 +
            m.f25*m.f25 + m.f34*m.f34 + m.f35*m.f35 + m.f45*m.f45)
m.v = Objective(rule=objective_function, sense=minimize)

m.lei_kf_1 = Constraint(expr= - m.f12 - m.f13 - m.f14 + m.g   ==  0)
m.lei_kf_2 = Constraint(expr=   m.f12 - m.f23 - m.f25         == 25)
m.lei_kf_3 = Constraint(expr=   m.f13 + m.f23 - m.f34 - m.f35 == 30)
m.lei_kf_4 = Constraint(expr=   m.f14 + m.f34 - m.f45         == 20)
m.lei_kf_5 = Constraint(expr=   m.f25 + m.f35 + m.f45         == 60)

m.f12_limit_a = Constraint(expr= m.f12 <=  100 * m.x12)
m.f12_limit_b = Constraint(expr= m.f12 >= -100 * m.x12)

m.f13_limit_a = Constraint(expr= m.f13 <=  100 * m.x13)
m.f13_limit_b = Constraint(expr= m.f13 >= -100 * m.x13)

m.f14_limit_a = Constraint(expr= m.f14 <=  100 * m.x14)
m.f14_limit_b = Constraint(expr= m.f14 >= -100 * m.x14)

m.f23_limit_a = Constraint(expr= m.f23 <=  100 * m.x23)
m.f23_limit_b = Constraint(expr= m.f23 >= -100 * m.x23)

m.f25_limit_a = Constraint(expr= m.f25 <=  100 * m.x25)
m.f25_limit_b = Constraint(expr= m.f25 >= -100 * m.x25)

m.f34_limit_a = Constraint(expr= m.f34 <=  100 * m.x34)
m.f34_limit_b = Constraint(expr= m.f34 >= -100 * m.x34)

m.f35_limit_a = Constraint(expr= m.f35 <=  100 * m.x35)
m.f35_limit_b = Constraint(expr= m.f35 >= -100 * m.x35)

m.f45_limit_a = Constraint(expr= m.f45 <=  100 * m.x45)
m.f45_limit_b = Constraint(expr= m.f45 >= -100 * m.x45)

m.g_limit     = Constraint(expr= inequality(0, m.g, 150))

m.condition_2 = Constraint(expr= m.x12 + m.x13 + m.x14 + m.x23 + m.x25 + m.x34 + m.x35 + m.x35 == 4)


solver = SolverFactory('cplex')
results = solver.solve(m, tee=False) 
print()
print('(Objective) v =  %d' % m.v())

print()

print('g = %.d' % m.g.value)

print()

print('x12 = %d' % m.x12.value)
print('x13 = %d' % m.x13.value)
print('x14 = %d' % m.x14.value)
print('x23 = %d' % m.x23.value)
print('x25 = %d' % m.x25.value)
print('x34 = %d' % m.x34.value)
print('x35 = %d' % m.x35.value)
print('x45 = %d' % m.x45.value)

print()

print('f12 = %d' % m.f12.value)
print('f13 = %d' % m.f13.value)
print('f14 = %d' % m.f14.value)
print('f23 = %d' % m.f23.value)
print('f25 = %d' % m.f25.value)
print('f34 = %d' % m.f34.value)
print('f35 = %d' % m.f35.value)
print('f45 = %d' % m.f45.value)
