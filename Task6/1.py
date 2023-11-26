# -*- coding: utf-8 -*-
import numpy as np
import math

def f(x):
    return math.sin(100 * x) * math.exp(-(x**2)) * math.cos(2 * x)

def main():
    h = 10**-5
    print("Метод трапеций: " + str(trap(f, 0, 3, h)))
    print("Метод Симпсона: " + str(simp(f, 0, 3, h)))

def trap(f, a, b, h):
    x = np.arange(a, b + h, h)
    su = 0
    for i in range(len(x) - 1):
        su += f(x[i]) + f(x[i + 1])
    return su * h / 2

def simp(f, a, b, h):
    x = np.arange(a, b + h, h)
    su = 0
    for i in range(len(x) - 1):
        su += f(x[i]) + f(x[i + 1]) + 4 * f(x[i] + h / 2)
    return su * h / 6

if __name__ == "__main__":
    main()
