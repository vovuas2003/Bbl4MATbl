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
    #A = np.asarray([[2.2,1,0.5,2],[1,1.3,2,1],[0.5,2,0.5,1.6],[2,1,1.6,2]])
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
    #c = np.zeros(n)
    #c[0] = 1
    #c = np.ones(n)
    c = np.asarray([n * (i + 1) for i in range(n)])
    B = np.zeros((n, n))
    for i in range(n):
        B[i, n - 1] = c[i]
    for i in range(n - 2, -1, -1):
        c = np.dot(A, c)
        for j in range(n):
            B[j, i] = c[j]
    f = np.dot(A, c)
    #l, _ = np.linalg.eig(A)
    #print(sorted(l))
    c = np.linalg.solve(B, f)
    c = [-1] + [x for x in c]
    #print(sorted(np.roots(c)))
    #########################################################
    print("\nЕсли n кругов пересекаются, то n собственных значений лежат в объединении этих кругов.\nВ данной задаче все 6 кругов пересекаются.\n")
    left = min([S[i][0] for i in range(n)])
    right = max([S[i][1] for i in range(n)])
    step = 10**-5
    big_step = 10**-4
    eps = 10**-4
    otv = []
    temp = left
    X = []
    Y = []
    y = 0
    while temp <= right:
        X.append(temp)
        y = -1 * polynom(c, temp)
        Y.append(y)
        if abs(y) <= eps:
            otv.append(temp)
            temp += big_step
        temp += step
    return 2 / (min(otv) + max(otv)), X, Y, otv

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

#Идея доработки Gersh_Kr(A) - вторая теорема Гершгорина.
#Если k кругов образуют связную область, изолированную от остальных n-k кругов,
#то первая область содержит ровно k, а вторая — n-k собственных значений матрицы A.
'''
Сохранять в подмассивах массива S ещё одно значение - количество решений на данном круге.
Сначала в каждом круге по одному решению.

Сортируем S по возрастанию координат левых концов кругов.

Проходимся по S, если правый конец текущего круга правее левого конца следующего,
то объединяем два этих круга в один (левый конец оставляем, правый - максимум из двух),
а количество решений складываем (к текущему прибавляем 1).
Если объединили круги, то рассматриваем объединённый круг и следующий.

Если нет пересечения, то переходим к следующей паре кругов.

В итоге если на круге одно решение, то можно находить на нём корень не тупым перебором,
а методом дихотомии (должно сработать; вроде если собственное значение имеет кратность,
то оно порождает пересекающиеся круги) или методом Ньютона.

Тупой перебор от левого конца круга до правого применять только если на круге несколько корней.
'''
