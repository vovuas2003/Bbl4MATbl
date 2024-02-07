import numpy as np
import matplotlib.pyplot as plt
import math

def f(t, u, eps):
    return 1 / (1 - u * math.cos(t)) + eps

def main():
    h = 10**-3
    l = -1.5
    s = 0
    r = 0.55
    t = np.arange(l, r, h)
    n = len(t)
    u = []
    plt.figure()
    u = Solve(t, n, h, f, 0, s)
    plt.plot(t, u)
    #
    e = np.arange(-0.005, 0.005, 0.0001)
    e = [0.001]
    e=[]
    #
    for eps in e:
        u = Solve(t, n, h, f, eps, s)
        '''
        if max(abs(u)) > 800:
            print(eps)
        '''
        plt.plot(t, u)
    plt.xlabel("t")
    plt.ylabel("u(t)")
    plt.tight_layout()
    plt.grid()
    plt.show()

def Solve(t, n, h, f, eps, s):
    u = [s]
    u.append(u[0] + h * f(t[0], u[0], eps))
    for i in range(n - 2):
        u.append(Next(u[-2], u[-1], t[i], t[i + 1], h, f, eps))
    u = np.asarray(u)
    return u

def Next(u_, u, t_, t, h, f, eps):
    return u_ + h * (f(t, u, eps) + f(t_, u_, eps))

if __name__ == "__main__":
    main()
