# -*- coding: utf-8 -*-
import numpy as np
import math

def f(x):
    return math.sin(x) / math.sqrt(x)

def f1(x):
    '''
    if np.isclose(x, 0):
        return math.sqrt(x)
    '''
    if (abs(x) <= 10**-4):
        return math.sqrt(x)
    return math.sin(x) / math.sqrt(x)

def f2(t):
    return 2 * math.sin(t**2)

def main():
    h = 10**-4
    zero = 2**-1
    eps = 10**-5
    i1 = trap(f, zero, 1, h)
    i2 = trap(f, zero / 2, 1, h)
    while(abs(i1 - i2) > eps):
        zero /= 2
        i1 = i2
        i2 = i2 = trap(f, zero / 2, 1, h)
    print("Метод трапеций, 3 способа устранения особенности:")
    print()
    print(str(trap(f, zero, 10, h)) + " (нижний предел -> 0)")
    print(str(trap(f1, 0, 10, h)) + " (sin(x) -> x при x -> 0)")
    print(str(trap(f2, 0, math.sqrt(10), h / 10)) + " (замена t = корень из x)")
    print()

def trap(f, a, b, h):
    x = np.arange(a, b + h, h)
    su = 0
    for i in range(len(x) - 1):
        su += f(x[i]) + f(x[i + 1])
    return su * h / 2

if __name__ == "__main__":
    main()
