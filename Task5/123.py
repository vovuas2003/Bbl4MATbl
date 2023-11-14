import matplotlib.pyplot as plt
import numpy as np

def main():
    t1 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    x1 = [1, 0.8, 0.5, 0.307, 0.2, 0.137, 0.1, 0.075, 0.06, 0.047, 0.039]
    newton_diff = []
    for i in range(1, len(t1)):
        newton_diff.append(diff(t1, x1, i))
    t1_ = np.linspace(t1[0], t1[-1], 10**3)
    x1_ = [New(t1, x1, newton_diff, t0) for t0 in t1_]
    plt.figure(figsize = (13.5, 6.3))
    plt.subplot(1, 3, 1)
    plt.plot(t1_, x1_)
    plt.title("Newton x(t)")
    plt.scatter(t1, x1)
    plt.grid()
    plt.tight_layout()

    t2 = [-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    y2 = [0.02, 0.079, 0.175, 0.303, 0.459, 0.638, 0.831, 1.03, 1.23, 1.42]
    t2_ = np.linspace(t2[0], t2[-1], 10**3)
    y2_ = [Lag(t2, y2, t0) for t0 in t2_]
    plt.subplot(1, 3, 2)
    plt.plot(t2_, y2_)
    plt.title("Lagrange y(t)")
    plt.scatter(t2, y2)
    plt.grid()
    plt.tight_layout()

    t = np.linspace(max(t1[0], t2[0]), min(t1[-1], t2[-1]), 10**3)
    x = [New(t1, x1, newton_diff, t0) for t0 in t]
    y = [Lag(t2, y2, t0) for t0 in t]
    xx = 0.431
    eps = 10**-3
    h = 10**-8
    ro = 3
    flag = False
    for t0 in t:
        if flag:
            if (abs(New(t1, x1, newton_diff, t0) - xx)) > h:
                tt2 = t0
                break
        if (abs(New(t1, x1, newton_diff, t0) - xx) < eps) and (flag == False):
            tt1 = t0
            flag = True
            continue
    ans = "dy/dx(" + str(xx) + ") = " + str(round((Lag(t2, y2, tt1) - Lag(t2, y2, tt2))/ (New(t1, x1, newton_diff, tt1) - New(t1, x1, newton_diff, tt2)), ro))
    plt.subplot(1, 3, 3)
    plt.plot(x, y)
    plt.title("y(x), " + ans)
    plt.grid()
    plt.tight_layout()
    plt.show()

def diff(x, y, k):
    r = 0
    for i in range(k + 1):
        d = 1
        for j in range(k + 1):
            if j != i:
                d *= x[i] - x[j]
        r += y[i] / d
    return r

def New(x, y, di, x0):
    r = y[0]
    for k in range(1, len(x)):
        p = 1
        for j in range(k):
            p *= (x0 - x[j])
        r += di[k - 1] * p
    return r

def lag_k(x, x0, k):
    p = 1
    q = 1
    xk = x[k]
    for i in range(len(x)):
        if i != k:
            p *= (x0 - x[i])
            q *= (xk - x[i])
    return p / q

def Lag(x, y, x0):
    r = 0
    for i in range(len(x)):
        r += lag_k(x, x0, i) * y[i]
    return r

if __name__ == "__main__":
    main()
