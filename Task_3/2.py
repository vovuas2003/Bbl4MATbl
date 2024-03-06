import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x**3

def main():
    N = 1000
    left = 0
    right = 1
    X = right - left
    h = X / N
    x = np.linspace(left, right, N)
    w = omega(X, N, x)
    g = [f(x[i]) for i in range(1, N)]
    c = np.linalg.solve(w, g)
    c = [i for i in c]
    l = lam(X, N)
    c = [c[i] / l[i] for i in range(N - 1)]
    y = [0]
    s = 0
    for i in range(N - 2):
        s = 0
        for k in range(N - 1):
            s -= c[k] * w[k][i] #pdf: s += ... ? 
        y.append(s)
    y.append(0)
    plt.figure(figsize = (13.5, 6.3))
    plt.suptitle("y'' = x^3, y(0) = y(1) = 0")
    plt.subplot(1, 2, 1)
    plt.title("Численное решение методом Фурье")
    plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.tight_layout()
    plt.grid()
    plt.subplot(1, 2, 2)
    plt.title("Точное аналитическое решение")
    #plt.plot(x, y, '--')
    plt.plot(x, [func(i) for i in x])
    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.tight_layout()
    plt.grid()
    plt.show()

def omega(X, N, x):
    w = []
    temp = []
    for n in range(1, N):
        temp = []
        for k in range(1, N):
            temp.append(((2 / X)**0.5) * np.sin(np.pi * k * x[n] / X))
        w.append(temp)
    return w

def lam(X, N):
    l = []
    for k in range(1, N):
        l.append((4 * ((np.sin((np.pi * k) / (2 * N)))**2)) / ((X / N)**2))
    return l

def func(x):
    return x * (x**4 - 1) / 20

if __name__ == "__main__":
    main()
