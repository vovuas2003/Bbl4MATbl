import numpy as np
import matplotlib.pyplot as plt

A = np.asarray([[0.78, 0.563], [0.457, 0.33]])
f = np.asarray([0.217, 0.127 + 0.0005])

def decLU(A):
    size = A.shape[0]
    U = np.zeros((size, size))
    L = np.zeros((size, size))
    for i in np.arange(0, size, 1):
        for j in np.arange(0, size, 1):
            if i == 0:
                U[i][j] = A[i][j]
                L[j][i] = A[j][i] / U[0][0]
            else:
                S = 0.0
                for k in np.arange(0, i - 1, 1):
                    S += L[i][k] * U[k][i]
                U[i][i] = A[i][i] - S
                S = 0.0
                for k in np.arange(0, i, 1):
                    S += L[i][k] * U[k][j]                    
                U[i][j] = A[i][j] - S
                S = 0.0
                for k in np.arange(0, i, 1):
                    S += L[j][k] * U[k][i]
                L[j][i] = (A[j][i] - S) / U[i][i]
    return U, L

def Solve(A, f):
    global decLU
    size = A.shape[0]
    y = np.zeros(size)
    x = np.zeros(size)
    U, L = decLU(A)
    for k in np.arange(0, size, 1):
        #S = 0.0
        #for j in np.arange(0, k-1, 1):
            #S += L[k][j] * u[j]
        y[k] = f[k] - np.dot(L[k][0 : k], y[0 : k])
    for k in np.arange(size - 1, -1, -1):
        x[k] = (y[k] - np.dot(U[k][k + 1 : size], x[k + 1 : size])) / U[k][k]
    return x

x = Solve(A, f)
print(x)
