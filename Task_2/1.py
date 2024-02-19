import numpy as np
import matplotlib.pyplot as plt
import math

'''
d2u / dt2 = f = -sin(u)
<=>
du / dt = y
dy / dt = f
'''

def f(t, x):
    u = x[0]
    du = x[1]
    
    func = -math.sin(u)
    
    return [x[1], func]

def main():
    x = [[1, 0]] #[[u0, y0] , [u1, y1], ...]
    h = 10**-3
    l = 0
    r = 4 * math.pi
    A = [[0.5 - math.sqrt(3) / 6, 0], [math.sqrt(3) / 3, -0.5 - math.sqrt(3) / 2]]
    B = [1 + math.sqrt(3) / 6, -math.sqrt(3) / 6]
    C = [0.5 - math.sqrt(3) / 6, -0.5 - math.sqrt(3) / 6]
    t = np.arange(l, r, h)
    n = len(t)
    for i in range(n - 1):
        x.append(Runge(t[i], x[-1], h, f, A, B, C))
    y = [x[i][0] for i in range(n)]
    a = []
    b = []
    left = -2.5
    right = 2.5
    step = 0.01
    c = np.arange(left, right + step, step)
    ii = 1j
    for i in c:
        for j in c:
            if i > step:
                continue
            if abs(R(i + j * ii, A, B)) <= 1:
                a.append(i)
                b.append(j)
    plt.figure(figsize = (13.5, 6.3))
    plt.suptitle("Метод абсолютно устойчив, но не А- и не L-устойчив")
    plt.subplot(1, 2, 1)
    plt.plot(t, y)
    plt.title("Решение уравнения")
    plt.xlabel("t")
    plt.ylabel("u(t)")
    plt.tight_layout()
    plt.grid()
    plt.subplot(1, 2, 2)
    plt.scatter(a, b)
    plt.title("Область абсолютной устойчивости")
    plt.xlim(left, right)
    plt.ylim(left, right)
    plt.xlabel("Re(z)")
    plt.ylabel("Im(z)")
    plt.tight_layout()
    plt.grid()
    plt.show()

def Runge(x, u, h, f, a, b, c):
    n = len(u)
    s = len(b)
    k = MPI(f, x, c, h, u, a, n, s)
    res = []
    for i in range(n):
        res.append(u[i] + h * sum([k[i][j] * b[j] for j in range(s)]))
    return res

def MPI(f, x, c, h, u, a, n, s):
    eps = 10**-3
    kold = []
    for i in range(n):
        temp = []
        for j in range(s):
            temp.append(0)
        kold.append(temp)
    k = []
    for i in range(n):
        temp = []
        for j in range(s):
            temp.append(0)
        k.append(temp)
    for i in range(s):
        U = [u[I] + h * sum([a[i][j] * kold[I][j] for j in range(s)]) for I in range(n)]
        F = f(x + c[i] * h, U)
        for j in range(n):
            k[j][i] = F[j]
    while(max([abs(k[i][j] - kold[i][j]) for i in range(n) for j in range(s)]) > eps):
        for i in range(n):
            for j in range(s):
                kold[i][j] = k[i][j]
        for i in range(s):
            U = [u[I] + h * sum([a[i][j] * kold[I][j] for j in range(s)]) for I in range(n)]
            F = f(x + c[i] * h, U)
            for j in range(n):
                k[j][i] = F[j]
    return k

def R(z, a, b):
    n = len(b)
    a = np.array(a).reshape(n, n)
    b = np.array(b).reshape(1, n)
    e = np.eye(n)
    x = np.linalg.det(e - z * a + z * (np.array([1, 1]).reshape(n, 1) * b))
    y = np.linalg.det(e - z * a)
    return x / y

if __name__ == "__main__":
    main()
