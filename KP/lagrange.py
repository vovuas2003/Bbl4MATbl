def koeff(x, x0, k):
    p = 1
    q = 1
    xk = x[k]
    for i in range(len(x)):
        if i != k:
            p *= (x0 - x[i])
            q *= (xk - x[i])
    return p / q

def func(x, y, x0):
    r = 0
    for i in range(len(x)):
        r += koeff(x, x0, i) * y[i]
    return r

def plot(x, y):
    import matplotlib.pyplot as plt
    import numpy as np
    x1 = np.linspace(x[0], x[-1], 10**3)
    y1 = [func(x, y, i) for i in x1]
    plt.plot(x1, y1)
    plt.scatter(x, y)
    plt.title("Lagrange f(x)")
    plt.grid()
    plt.show()

def cheb(a, b, n):
    x = []
    import math
    for k in range(1, n + 1):
        x.append((a + b) / 2 + ((b - a) * math.cos(((2 * k - 1) * math.pi) / (2 * n))) / 2)
    return x
