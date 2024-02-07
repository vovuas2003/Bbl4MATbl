import numpy as np
import matplotlib.pyplot as plt
import math

def f(t, u, eps):
    return 1 / (1 - u * math.cos(t)) + eps

def main():
    h = 10**-3
    l = -1.5
    r = 0.55
    t = np.arange(l, r, h)
    n = len(t)
    u = Solve(t, n, h, f, 0)
    plt.figure()
    plt.plot(t, u)
    u = Solve(t, n, h, f, 10**-3)
    plt.plot(t, u)
    u = Solve(t, n, h, f, -10**-3)
    plt.plot(t, u)
    plt.tight_layout()
    plt.grid()
    plt.show()

def Solve(t, n, h, f, eps):
    u = [0, 0]
    for i in range(n - 2):
        u.append(Next(u[-2], u[-1], t[i], t[i + 1], h, f, eps))
    u = np.asarray(u)
    return u

def Next(u_, u, t_, t, h, f, eps):
    return u_ + h * (f(t, u, eps) + f(t_, u_, eps))

if __name__ == "__main__":
    main()
