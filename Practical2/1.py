# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math

def f(x):
    return 1 / (1 + 25 * x * x)

def main():
    a = -1
    b = 1
    x = np.linspace(a, b, 10**3)
    y = [f(i) for i in x]
    plt.figure(figsize = (13.5, 6.3))
    plt.subplot(1, 3, 1)
    plt.title("Исходный график")
    plt.plot(x, y)
    plt.tight_layout()
    plt.grid()
    N = [4, 6, 10]
    plt.subplot(1, 3, 2)
    plt.title("Узлы равномерно")
    plt.tight_layout()
    plt.grid()
    for n in N:
        x1 = np.linspace(a, b, n)
        y1 = [f(i) for i in x1]
        y = [Lagr(x1, y1, i) for i in x]
        plt.plot(x, y)
        plt.scatter(x1, y1)
    plt.subplot(1, 3, 3)
    plt.title("Узлы - нули полинома Чебышева")
    plt.tight_layout()
    plt.grid()
    for n in N:
        x1 = cheb(a, b, n)
        y1 = [f(i) for i in x1]
        F = diff(x1, y1)
        y = [Newt(x1, F, i) for i in x]
        plt.plot(x, y)
        plt.scatter(x1, y1)
    plt.show()

def koeff(x, x0, k):
    p = 1
    q = 1
    xk = x[k]
    for i in range(len(x)):
        if i != k:
            p *= (x0 - x[i])
            q *= (xk - x[i])
    return p / q

def Lagr(x, y, x0):
    r = 0
    for i in range(len(x)):
        r += koeff(x, x0, i) * y[i]
    return r

def cheb(a, b, n):
    x = []
    for k in range(1, n + 1):
        x.append((a + b) / 2 + ((b - a) * math.cos(((2 * k - 1) * math.pi) / (2 * n))) / 2)
    return x

def diff(x, y):
    n = len(x)
    F = [y]
    df = 0
    dx = 0
    temp = []
    for i in range(1, n):
        temp = []
        for j in range(n - i):
            df = F[i - 1][j + 1] - F[i - 1][j]
            dx = x[j + i] - x[j]
            temp.append(df / dx)
        F.append(temp)
    return F

def Newt(x, F, x0):
    s = 0
    p = 1
    for i in range(len(F)):
        p = 1
        for j in range(i):
            p *= (x0 - x[j])
        s += F[i][0] * p
    return s

if __name__ == "__main__":
    main()
