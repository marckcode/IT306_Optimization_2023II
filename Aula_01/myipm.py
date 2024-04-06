import numpy as np
from numpy.linalg import *

def myipm(c, A, b):
    # passo 1
    m, n = A.shape
    erro = 10**-5
    gamma = 0.99995
    
    # passo 2
    xtil = np.zeros((n, 1))
    
    for j in range(0, n):
        xtil[j] = 1 / (norm(A[:, j], 2) + 1)
        
    eta = (norm(b, 2) + 1) / (norm(A.dot(xtil), 2) + 1)
    
    x = eta * xtil
    y = np.zeros((m, 1))
    z = np.zeros((n, 1))
    
    print("x = ")
    print(x)
    print()
    
    print("y = ")
    print(y)
    print()
    
    for j in range(0, n):
        if c[j] < 1:
            z[j] = 1
        else:
            z[j] = c[j]
            
    print("z = ")
    print(z)
    print()
    
    mu = np.dot(z.T, x) / n
    print("mu = %.4f" % mu[0][0])
    
    input() # Pressione "enter" para continuar
    
    # passo 3
    fp = norm(A.dot(x)[:,0] - b, 2) / (1 + norm(x, 2))
    fd = norm(A.T.dot(y)[:, 0] + z[:,0] - c, 2) / (1 + norm(y, 2) + norm(z, 2))
    fo = norm(c.T.dot(x) - b.T.dot(y), 2) / (1 + norm(b.T.dot(y), 2))
    print("fp = %.4f" % fp)
    print("fd = %.4f" % fd)
    print("fo = %.4f" % fo)
    
    input() # Pressione "enter" para continuar
    itr = 0
    while fp > erro or fd > erro or fo > erro:
        itr = itr + 1
        print(f"============== ITER #{itr}")
        
        # passo 4
        Jnewton = np.block([[                 A,  np.zeros((m, m)),     np.zeros((m, n))  ],
                            [  np.zeros((n, n)),               A.T,            np.eye(n)  ],
                            [   np.diag(z[:,0]),  np.zeros((n, m)),      np.diag(x[:,0])  ]])

        print("Jnewton = ")
        print(Jnewton)
        print()
        prod = np.diag(x[:,0]).dot(np.diag(z[:,0])).dot(np.ones((n,1)))
    
        m1 = (b - A.dot(x[:,0]))
        m2 = (c - A.T.dot(y[:,0]) - z[:,0])
        m3 = (mu * np.ones((n, 1)) - prod).squeeze()

        bnewton = np.block([[m1, m2, m3]]).T
        print("bnewton = ")
        print(bnewton)
        
        input() # Pressione "enter" para continuar
        
        # passo 5
        dw = solve(Jnewton, bnewton)
        print("dw = ")
        print(dw)

        dx = dw[:n]
        dy = dw[n:n+m]
        dz = dw[n+m:2*n+m]
        
        input() # Pressione "enter" para continuar
        
        # passo 6
        alphap = 1
        for i in range(0, n):
            if dx[i][0] < 0:
                if -x[i][0]/dx[i][0] < alphap:
                    alphap = -x[i][0]/dx[i][0]
        
        print("alphap = ", alphap)
        print()
        
        alphad = 1
        for i in range(0, n):
            if dz[i][0] < 0:
                if -z[i][0]/dz[i][0] < alphad:
                    alphad = -z[i][0]/dz[i][0]
        
        print("alphad = ", alphad)

        input() # Pressione "enter" para continuar
        
        # passo 7
        x = x + gamma * alphap * dx
        y = y + gamma * alphad * dy
        z = z + gamma * alphad * dz
        
        print("x = ")
        print(x)
        print()
        
        print("y = ")
        print(y)
        print()
        
        print("z = ")
        print(z)
        
        input() # Pressione "enter" para continuar
        
        # passo 8
        if c.T.dot(x)[0] >= b.T.dot(y)[0]:
            sigma = 0.1
        else:
            sigma = 2
        
        print("sigma = %.4f" % sigma)
        print()
        
        mu = (sigma * (z.T.dot(x)[0])/n)[0]
        if mu < 10**-5:
            mu = 10**-5
        print("mu = %.4f" % mu)
        
        input() # Pressione "enter" para continuar
    
        # passo 3
        fp = norm(A.dot(x)[:,0] - b, 2) / (1 + norm(x, 2))
        fd = norm(A.T.dot(y)[:, 0] + z[:,0] - c, 2) / (1 + norm(y, 2) + norm(z, 2))
        fo = norm(c.T.dot(x) - b.T.dot(y), 2) / (1 + norm(b.T.dot(y), 2))
        print("fp = ", fp)
        print("fd = ", fd)
        print("fo = ", fo)

        input() # Pressione "enter" para continuar
    
    print("# de iteraçoes: %d" % itr)
    print("ans = ")
    print(x)

# Definir os dados
c = np.array([-2, -3, 0, 0])
A = np.array([[2, 1, 1, 0], [1, 2, 0, 1]])
b = np.array([8, 6])

myipm(c, A, b)