import numpy as np
import matplotlib.pyplot as plt
import math

def f(u):
    return -math.sin(u)

def main():
    h = 10**-3
    l = 0
    r = 4 * math.pi
    t = np.arange(l, r, h)
    n = len(t)
    u1 = [1]
    for i in range(n - 1):
        u1.append(Euler(u1[-1], h, f))
    plt.figure()
    plt.plot(t, u1)
    #plt.plot(t, u2)
    #plt.legend
    plt.grid()
    plt.tight_layout()
    plt.show()

def Euler(u, h, f):
    return u + h * f(u)

if __name__ == "__main__":
    main()
