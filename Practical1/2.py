# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import math

def main():
    n = 6
    A = np.zeros((n, n))
    tay = 1
    eps = 10**-6
    for i in range(n):
        A[i, i] = 2 + (1 / (n**2))
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
    if printxr:
        print("xr = ", end = '')
        print(xr)
        print()
    x, y = solve(A, f, tay, eps)
    n = norm(x - xr)
    if printx:
        print("x = ", end = '')
        print(x)
        print()
    if printn:
        print("norm = ", end = '')
        print(n)
        print()
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
    tay = Gersh_Kr(A)
    x, y = solve(A, f, tay, eps)
    n = norm(x - xr)
    if printx:
        print("x = ", end = '')
        print(x)
        print()
    if printn:
        print("norm = ", end = '')
        print(n)
        print()
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
    if printx:
        print("x = ", end = '')
        print(x)
        print()
    if printn:
        print("norm = ", end = '')
        print(n)
        print()
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
        S.append([A[i, i], r])
    return 1

def by_numpy(A):
    l, _ = np.linalg.eig(A)
    t = 2 / (max(l) + min(l))
    return t

def norm(x):
    return math.sqrt(sum([abs(i)**2 for i in x]))

if __name__ == "__main__":
    main()
