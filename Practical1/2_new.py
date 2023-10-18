# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import math

def main():
    n = 6
    A = np.zeros((n, n))
    tay = 0.4
    eps = 10**-4
    for i in range(n):
        A[i, i] = 2 + (i / n)**2
    for i in range(1, n):
        A[i, i - 1] = -1
    for i in range(n - 1):
        A[i, i + 1] = -1
    A[0, n - 1] = -1
    A[n - 1, 0] = -1
    f = np.zeros(n)
    for i in range(n):
        f[i] = (1 + (n**2)*((math.sin(math.pi / n))**2)) * math.sin((2 * math.pi * (i)) / n)
    global_solve(A, f, tay, eps)

def global_solve(A, f, tay, eps):
    zoom = False
    textx = False
    printx = True
    textn = False
    printn = True
    printxr = True
    xc = 4.2
    yc = 0.37
    xn = 2.5
    yn = 0.45
    xr = np.linalg.solve(A, f)
    Printxr(printxr, xr)
    x, y = solve(A, f, tay, eps)
    n = norm(x - xr)
    Printx(printx, x)
    Printn(printn, n)
    plt.figure(figsize = (13.5, 6.3))
    plt.subplot(1, 3, 1)
    plt.plot([int(i + 1) for i in range(len(y))], y)
    if textx:
        plt.text(xc, yc, "x = \n\n" + "\n".join([str(round(i, 2)) for i in x]))
    if textn:
        plt.text(xn, yn, "norm = " + str(round(n, 3)))
    plt.title("Введённый вручную параметр")
    plt.xlabel("Номер итерации")
    plt.ylabel("Третья норма невязки Ax - f")
    plt.tight_layout()
    tay, X, Y, otv = Gersh_Kr(A)
    x, y = solve(A, f, tay, eps)
    n = norm(x - xr)
    Printx(printx, x)
    Printn(printn, n)
    plt.subplot(1, 3, 2)
    plt.plot([int(i + 1) for i in range(len(y))], y)
    if textx:
        plt.text(xc, yc, "x = \n\n" + "\n".join([str(round(i, 2)) for i in x]))
    if textn:
        plt.text(xn, yn, "norm = " + str(round(n, 3)))
    plt.title("Параметр методом Крылова из кругов Гершгорина")
    plt.xlabel("Номер итерации")
    plt.ylabel("Третья норма невязки Ax - f")
    plt.tight_layout()
    tay = by_numpy(A)
    x, y = solve(A, f, tay, eps)
    n = norm(x - xr)
    Printx(printx, x)
    Printn(printn, n)
    plt.subplot(1, 3, 3)
    plt.plot([int(i + 1) for i in range(len(y))], y)
    if textx:
        plt.text(xc, yc, "x = \n\n" + "\n".join([str(round(i, 2)) for i in x]))
    if textn:
        plt.text(xn, yn, "norm = " + str(round(n, 3)))
    plt.title("Параметр с помощью np.linalg.eig")
    plt.xlabel("Номер итерации")
    plt.ylabel("Третья норма невязки Ax - f")
    plt.tight_layout()
    if zoom:
        try:
            plt.get_current_fig_manager().window.state('zoomed')
        except:
            print("\nОшибка при попытке развернуть графики на весь экран!")
    plt.show()
    plt.figure()
    plt.plot(X, Y)
    plt.title("Количество корней характеристического многочлена: " + str(len(otv)))
    plt.scatter(otv, [0] * len(otv))
    plt.grid()
    plt.show()

def solve(A, f, tay, eps):
    max_i = 10**5
    n = len(A)
    x = np.zeros(n)
    y = []
    nev = 0
    E = np.eye(n)
    for i in range(max_i):
        x = (E - tay * A).dot(x) + tay * f
        nev = norm(A.dot(x) - f)
        y.append(nev)
        if nev < eps:
            break
    return x, y

def Gersh_Kr(A):
    n = len(A)
    S = []
    r = 0
    for i in range(n):
        r = 0
        for j in range(n):
            if i == j:
                continue
            r += abs(A[i, j])
        S.append([A[i, i] - r, A[i, i] + r])
    print("Оценка собственных значений с помощью кругов Гершгорина:")
    for i in range(n):
        print("lambda " + str(i + 1) + " is in [" + str(S[i][0]) + "; " + str(S[i][1]) + "]")
    c = np.asarray([n * (i + 1) for i in range(n)])
    B = np.zeros((n, n))
    for i in range(n):
        B[i, n - 1] = c[i]
    for i in range(n - 2, -1, -1):
        c = np.dot(A, c)
        for j in range(n):
            B[j, i] = c[j]
    f = np.dot(A, c)
    c = np.linalg.solve(B, f)
    c = [-1] + [x for x in c]
    print("\nЕсли n кругов пересекаются, то n собственных значений лежат в объединении этих кругов.\nВ данной задаче все 6 кругов пересекаются.\n")
    left = min([S[i][0] for i in range(n)])
    right = max([S[i][1] for i in range(n)])
    otv = []
    X = np.linspace(left, right, 10**5)
    Y = [-1 * polynom(c, i) for i in X]
    step = (right - left) / (5 * n)
    x_prev = left
    x_temp = x_prev + step
    while x_temp <= right:
        if polynom(c, x_prev) * polynom(c, x_temp) < 0:
            otv.append(dihot(c, x_prev, x_temp))
        x_prev = x_temp
        x_temp += step
    return 2 / (min(otv) + max(otv)), X, Y, otv

def dihot(c, x_prev, x_temp):
    eps = 10**-3
    x = (x_temp + x_prev) / 2
    while (x_temp - x_prev) >= eps: 
        x = (x_temp + x_prev) / 2
        if polynom(c, x_prev) * polynom(c, x) < 0:
            x_temp = x
        else:
            x_prev = x
    return (x_temp + x_prev) / 2

def polynom(c, x):
    p = 0
    for i in c:
        p = p * x + i
    return p

def by_numpy(A):
    l, _ = np.linalg.eig(A)
    t = 2 / (max(l) + min(l))
    return t

def norm(x):
    return math.sqrt(sum([abs(i)**2 for i in x]))

def Printxr(f, xr):
    if f:
        print("xr = ", end = '')
        print(xr)
        print()

def Printx(f, x):
    if f:
        print("x = ", end = '')
        print(x)
        print()

def Printn(f, n):
    if f:
        print("norm = ", end = '')
        print(n)
        print()

if __name__ == "__main__":
    main()
