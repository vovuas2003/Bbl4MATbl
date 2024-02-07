import numpy as np
import matplotlib.pyplot as plt
import math

#функция не зависит от t, поэтому методы Эйлера и Рунге-Кутты
#записаны для этого частного случая
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
    u1 = np.asarray(u1)
    u2 = [1]
    for i in range(n - 1):
        u2.append(RungeKutt(u2[-1], h, f))
    u2 = np.asarray(u2)
    plt.figure(figsize = (13.5, 6.3))
    plt.suptitle("max|u1(t) - u2(t)| = " + str(round(max(abs(u1 - u2)), 6)))
    plt.subplot(1, 2, 1)
    plt.title("Метод Эйлера u1(t)")
    plt.xlabel("t")
    plt.ylabel("u1(t)")
    plt.plot(t, u1)
    plt.tight_layout()
    plt.grid()
    plt.subplot(1, 2, 2)
    plt.title("Метод Рунге-Кутты u2(t)")
    plt.xlabel("t")
    plt.ylabel("u2(t)")
    plt.plot(t, u2)
    plt.tight_layout()
    plt.grid()
    plt.show()

def Euler(u, h, f):
    return u + h * f(u)

def RungeKutt(u, h, f):
    k1 = f(u)
    k2 = f(u + h * k1 / 2)
    k3 = f(u + h * k2 / 2)
    k4 = f(u + h * k3)
    return u + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6

if __name__ == "__main__":
    main()
