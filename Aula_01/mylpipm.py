import numpy as np
from numpy.linalg import *

def mylpipm(A, b, c):
    # passo 1
    ef = 10**-5
    eo = 10**-5
    gamma = 0.99995
    
    # passo 2
    k = 0
    m, n = A.shape
    xch = np.zeros((n, 1))
    
    for j in range(0, n):
        xch[j] = 1 / (norm(A[:, j], 2) + 1)

    nu = (norm(b, 2) + 1) / (norm(A.dot(xch), 2) + 1)
    
    xk = nu * xch
    zk = np.zeros((n, 1))
    
    for j in range(0, n):
        if c[j] < 1:
            zk[j] = 1
        else:
            zk[j] = c[j]
            
    yk = np.zeros((m, 1))

    muk = zk.T.dot(xk[:,0])[0] / n
    
    # passo 3
    fp = norm(A.dot(xk)[:,0] - b, 2) / (1 + norm(xk, 2))
    fd = norm(A.T.dot(yk)[:, 0] + zk[:,0] - c, 2) / (1 + norm(yk, 2) + norm(zk, 2))
    fo = norm(c.T.dot(xk) - b.T.dot(yk), 2) / (1 + norm(b.T.dot(yk), 2))
    
    itr = 0
    while fp > ef or fd > ef or fo > eo:
        itr = itr + 1
        
        # passo 4
        MJ = np.block([[                   A, np.zeros((m, m)),     np.zeros((m, n))],
                            [    np.zeros((n, n)),              A.T,            np.eye(n)],
                            [     np.diag(zk[:,0]), np.zeros((n, m)),      np.diag(xk[:,0])]])

        prod = np.diag(xk[:,0]).dot(np.diag(zk[:,0])).dot(np.ones((n,1)))
        m1 = (b - A.dot(xk[:,0]))
        m2 = (c - A.T.dot(yk[:,0]) - zk[:,0])
        m3 = (muk*np.ones((n,1))- prod).squeeze()

        Mb = np.block([[m1, m2, m3]]).T

        # passo 5
        dw = solve(MJ, Mb)

        dx = dw[:n]
        dy = dw[n:n+m]
        dz = dw[n+m:2*n+m]
        
        # passo 6
        alphap = 1
        for i in range(0, n):
            if dx[i][0] < 0:
                alphapaux = -xk[i][0]/dx[i][0]
                if alphapaux < alphap:
                    alphap = alphapaux
            
        alphad = 1
        for i in range(0, n):
            if dz[i][0] < 0:
                alphadaux = -zk[i][0]/dz[i][0]
                if alphadaux < alphad:
                    alphad = alphadaux

        # passo 7
        xk = xk + gamma * alphap * dx
        yk = yk + gamma * alphad * dy
        zk = zk + gamma * alphad * dz
        
        # passo 8
        if c.T.dot(xk)[0] >= b.T.dot(yk)[0]:
            sigma = 0.1
        else:
            sigma = 2
            
        muk = (sigma * (zk.T.dot(xk)[0])/n)[0]
        if muk < 10**-5:
            muk = 10**-5
            
        # passo 9
        k = k + 1
        
        fp = norm(A.dot(xk)[:,0] - b, 2) / (1 + norm(xk, 2))
        fd = norm(A.T.dot(yk)[:, 0] + zk[:,0] - c, 2) / (1 + norm(yk, 2) + norm(zk, 2))
        fo = norm(c.T.dot(xk) - b.T.dot(yk), 2) / (1 + norm(b.T.dot(yk), 2))

    x = xk
    fo = c.T.dot(x[:, 0])
    
    print("# de iteraçoes = ", itr)
    print("ans = ")
    print(x)
    print("fo = ", fo)
    print()

# Definir os dados
A = np.array([[2, 1, 1, 0], [1, 2, 0, 1]])
b = np.array([8, 6])
c = np.array([-2, -3, 0, 0])
mylpipm(A, b, c)