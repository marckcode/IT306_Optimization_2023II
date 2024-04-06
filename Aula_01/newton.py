#import math
from math import *

# passo 1
k = 0
print("k = %d" % k)
erro = 0.001
xk = 1
print("xk = %d" % xk)

# passo 2
r = xk + 5 * cos(xk)
print("r = %.4f" % r)

# passo 3
while abs(r) > erro:
    # passo 4
    dr = 1 - 5 * sin(xk)
    print("dr = %.4f" % dr)
    
    # passo 5
    dx = - r / dr
    print("dx = %.4f" % dx)
    xk = xk + dx
    print("xk = %.4f" % xk)
    
    # passo 6
    k = k + 1
    print("\nk = %d" % k)
    
    # passo 2
    r = xk + 5 * cos(xk)
    print("r = %.4f" % r)


