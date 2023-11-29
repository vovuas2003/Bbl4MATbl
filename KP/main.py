import newton as new
import lagrange as lag
import integral as inte

import math
import matplotlib.pyplot as plt
import numpy as np

def n1():
    x = [-2.5, -2.17, -1.5, -1.17, -0.83, 2.5]
    y = [-0.52, -0.71, -0.77, -0.61, -0.31, 11.13]
    print(lag.func(y, x, 0))

def n2():
    x = [1.3, 2.4, 4.7, 6.2, 7.1]
    y = [0.2, -0.12, 0.45, 1.78, 1.34]
    new.table(x, y, 5)
    h = 10**-6
    n = 2
    x0 = 3.45
    print(new.deriv(x, y, h, n, x0))

def n3():
    def f(x):
        return 1 + math.cos(x) / (10 + x)
    a = -3
    b = 4
    n = 4
    x = lag.cheb(a, b, n)
    y = [f(i) for i in x]
    x1 = np.linspace(a, b, 10**3)
    y1 = [abs(lag.func(x, y, i) - f(i)) for i in x1]
    print(max(y1))
    lag.plot(x, y)

def n4():
    def f(x):
        return math.cos(x) / (2 + x * x)
    a = 0
    b = 10**4
    eps = 10**-3
    inte.trap_auto(f, a, b, eps)
    print(inte.trap(f, a, b, 1))

def n5():
    def f(x):
        return x**3 + 5 * (x**2) + 10 * x + 4.5
    a = 0
    b = 4
    print(inte.gauss3(f, a, b))

def main():
    n1()
    print()
    n2()
    print()
    n3()
    print()
    n4()
    print()
    n5()

if __name__ == "__main__":
    main()
