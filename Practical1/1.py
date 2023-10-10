# -*- coding: utf-8 -*-
import numpy as np
import math

def main():
    n = 6
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = 10 / ((i + 1) + (j + 1) - 1)
    f = np.zeros(n)
    for i in range(n):
        f[i] = sum([A[i, j] for j in range(n)])
    global_solve(A, f)

def global_solve(A, f):
    L = LL(A)
    print("Матрица L:")
    print(L)
    print()
    x = solve(L, f)
    print("Решение методом Холецкого:")
    print(x)
    print()
    x_np = np.linalg.solve(A, f)
    print("Решение с помощью np.linalg.solve:")
    print(x_np)
    print()
    print("Третья (Евклидова) норма разницы решений:")
    print(norm(x - x_np))

def LL(A):
    n = len(A)
    L = np.zeros((n, n))
    for i in range(n):
        for j in range(i):
            L[i, j] = (A[i, j] - sum([L[i, k] * L[j, k] for k in range(j)])) / L[j, j]
        L[i, i] = math.sqrt(A[i, i] - sum([L[i, k]**2 for k in range(i)]))
    return L

def solve(L, f):
    n = len(L)
    v = np.zeros(n)
    for i in range(n):
        v[i] = (f[i] - sum([L[i, k] * v[k] for k in range(i)])) / L[i, i]
    u = np.zeros(n)
    for i in range(n - 1, -1, -1):
        u[i] = (v[i] - sum([L[j, i] * u[j] for j in range(i, n)])) / L[i, i]
    return u

def norm(x):
    return math.sqrt(sum([abs(i)**2 for i in x]))

if __name__ == "__main__":
    main()
