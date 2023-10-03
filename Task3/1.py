# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math

#настройка
N = 30

#1-я норма матрицы
def n1(a):
    a = abs(a)
    su = a.sum(axis = 1) #массив сумм строк
    return max(su)
#2-я норма матрицы
def n2(a):
    a = abs(a)
    su = a.sum(axis = 0) #массив сумм столбцов
    return max(su)
#3-я норма матрицы
def n3(a):
    aa = np.dot(np.transpose(np.conj(a)), a)
    l, _ = np.linalg.eig(aa) #массив собственных чисел
    l = [abs(i) for i in l]
    return math.sqrt(max(l))
#создание матрицы
def make_A(n):
    A = np.eye(n)
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            A[i][j] = -1
    for j in range(n):
        A[n - 1][j] = 1
    return A
#LU разложение
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
#решение после LU разложения
def Solve(A, f):
    global decLU
    size = A.shape[0]
    y = np.zeros(size)
    x = np.zeros(size)
    U, L = decLU(A)
    for k in np.arange(0, size, 1):
        y[k] = f[k] - np.dot(L[k][0 : k], y[0 : k])
    for k in np.arange(size - 1, -1, -1):
        x[k] = (y[k] - np.dot(U[k][k + 1 : size], x[k + 1 : size])) / U[k][k]
    return x

x = []
y1 = []
y2 = []
y3 = []
for n in range(1, N + 1):
    f = [1] * n
    A = make_A(n)
    X = Solve(A, f)
    print("n = " + str(n))
    print(X)
    print()
    A_1 = np.linalg.inv(A)
    x.append(n)
    y1.append(n1(A) * n1(A_1))
    y2.append(n2(A) * n2(A_1))
    y3.append(n3(A) * n3(A_1))

plt.figure()
plt.subplot(1, 3, 1)
plt.title("1-я норма")
plt.plot(x, y1)
plt.xlabel("n")
plt.ylabel("Число обусловенности")
plt.tight_layout()
plt.subplot(1, 3, 2)
plt.title("2-я норма")
plt.plot(x, y2)
plt.xlabel("n")
plt.ylabel("Число обусловенности")
plt.tight_layout()
plt.subplot(1, 3, 3)
plt.title("3-я норма")
plt.plot(x, y3)
plt.xlabel("n")
plt.ylabel("Число обусловенности")
plt.tight_layout()
plt.show()
