# -*- coding: utf-8 -*-
import numpy as np
import math

def f(x):
    return math.cos(x) / (2 + x**2)

def main():
    eps = 10**-4
    h = 1
    infi = 10**5
    i1 = simp(f, 0, infi, h)
    i2 = simp(f, 0, infi, h / 2)
    while(abs(i1 - i2) > (eps * (2**4 - 1))):
        h /= 2
        i1 = i2
        i2 = simp(f, 0, infi, h / 2)
    print("I = " + str(i1))
    print("h = " + str(h))

def simp(f, a, b, h):
    x = np.arange(a, b + h, h)
    su = 0
    for i in range(len(x) - 1):
        su += f(x[i]) + f(x[i + 1]) + 4 * f(x[i] + h / 2)
    return su * h / 6

if __name__ == "__main__":
    main()
