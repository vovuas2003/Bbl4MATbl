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

def table(x, y, r = -1):
    F = diff(x, y)
    print("=" * 68)
    print("x", end = ' ')
    print(*x)
    print("f(x) == f0", end = ' ')
    print(*y)
    if r == -1:
        for i in range(1, len(F)):
            print("f" + str(i) + " ")
            for j in F[i]:
                print(j, end = ' ')
            print()
    else:
        for i in range(1, len(F)):
            print("f" + str(i) + " ")
            for j in F[i]:
                print(round(j, r), end = ' ')
            print()
    print("=" * 68)

def func(x, y, x0):
    F = diff(x, y)
    s = 0
    p = 1
    for i in range(len(F)):
        p = 1
        for j in range(i):
            p *= (x0 - x[j])
        s += F[i][0] * p
    return s

def func_diff(x, F, x0):
    s = 0
    p = 1
    for i in range(len(F)):
        p = 1
        for j in range(i):
            p *= (x0 - x[j])
        s += F[i][0] * p
    return s

def plot(x, y):
    import matplotlib.pyplot as plt
    import numpy as np
    F = diff(x, y)
    x1 = np.linspace(x[0], x[-1], 10**3)
    y1 = [func_diff(x, F, i) for i in x1]
    plt.plot(x1, y1)
    plt.scatter(x, y)
    plt.title("Newton f(x)")
    plt.grid()
    plt.show()

def fac(x):
    p = 1
    while x != 0:
        p *= x
        x -=1
    return p

def soch(n, k):
    return fac(n) // (fac(n - k) * fac(k))

def deriv_f(f, x, h, n):
    s = 0
    for k in range(0, n + 1):
        s += soch(n, k) * f(x + k * h) * ((-1) ** (n - k))
    return s / (h ** n)

def deriv(x, y, h, n, x0):
    F = diff(x, y)
    s = 0
    for k in range(0, n + 1):
        s += soch(n, k) * func_diff(x, F, x0 + k * h) * ((-1) ** (n - k))
    return s / (h ** n)
