# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math

A = np.asarray([[0.78, 0.563], [0.457, 0.33]])
f = np.asarray([0.217, 0.127 + 0.0005])
d1 = np.asarray([0, 0.0005])
d2 = np.asarray([0.0001, 0])
d3 = np.asarray([0.001, 0.0006])
x = np.asarray([1, -1])

#1-я норма матрицы
def norm_m(a):
    a = abs(a)
    su = a.sum(axis = 1) #массив сумм строк
    return max(su)
#число обусловленности
def mu(a):
    return norm_m(a) * norm_m(np.linalg.inv(a))
#1-я норма вектора
def norm_v(x):
    return max(abs(x))

def solve(A, f):
    lambd, _ = np.linalg.eig(A)
    t = 2 / (max(lambd) + min(lambd))
    #prev = np.asarray([A[i][i]/f[i] for i in range(len(f))])
    prev = np.zeros(len(f))
    x = (np.eye(len(f)) - t * A).dot(prev) + t * f
    while(norm_v(x - prev) > 10**-6):
        prev = x
        x = (np.eye(len(f)) - t * A).dot(prev) + t * f
    return x

x1 = solve(A, f + d1)
x2 = solve(A, f + d2)
x3 = solve(A, f + d3)

m = mu(A)
print(x1, norm_v(x - x1), (norm_v(x - x1) / norm_v(x)) <= (m * norm_v(d1) / norm_v(f + d1)))
print(x2, norm_v(x - x2), (norm_v(x - x2) / norm_v(x)) <= (m * norm_v(d2) / norm_v(f + d2)))
print(x3, norm_v(x - x3), (norm_v(x - x3) / norm_v(x)) <= (m * norm_v(d3) / norm_v(f + d3)))
