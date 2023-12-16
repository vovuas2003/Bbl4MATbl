#обоснование формул(стр. 24):
#https://kpfu.ru/portal/docs/F597883768/2017.Chislennye.metody.resheniya.IU.i.kompleks.programm.na.yazyke.Matlab.pdf
import numpy as np
import math
import Lagrange_method as LAGR
import Integral
import matplotlib.pyplot as plt
class I_FUNC:
    def __init__(self, a, b):
        self.a = a
        self.b = b
    def g(self):
        return 1
    def K(self, x_i, x_j):
        return 0.2/(0.04 + (x_i - x_j)**2)
    def f(self, x):
        return math.cos(math.pi * x)
    def l(self):
        return 1
    def a_b(self):
        return self.a ,self.b

def solve(F, n):
    # задаем параметры для вычислений матрицы А
    a, b = F.a_b()
    g = F.g()
    l = F.l()
    koef = l/g # нужен для домножения значений в А
    f = F.f
    F_K = F.K

    #задаем узлы для x и s
    nodes_array = list(np.arange(a, b, (b-a)/(n - 1))) + [b]

    # задаем матрицу решений f_values
    f_values = []
    for i in nodes_array:
        f_values.append(f(i)/koef)
    f_values = np.array(f_values)

    # Задаем квадратную матрицу K(x_i, x_j)
    K = []
    for i in nodes_array:
        k = []
        for j in nodes_array:
            k.append(F_K(i, j))
        K.append(k)
    K = np.array(K)

    # считаем таблицу весов
    A_I_J = []
    for i in range(len(nodes_array)):
        arr = []
        pol, basic_polinom = LAGR.create_Lagrange_polynomial(K[i], [])
        for j in range(len(basic_polinom)):
            arr.append(Integral.Simpson_method(-1, 1, basic_polinom[i]))
        A_I_J.append(arr)
    A_I_J = np.array(A_I_J)

    #составляем матрицу коэффициентов уравнения u(x)
    SOLVE = []
    for i in range(len(nodes_array)):
        solve = []
        for j in range(len(nodes_array)):
            if j == i:
                solve.append(1 - ((b - a) / 2 ) * A_I_J[i][j] * K[i][j] / koef)
            else:
                solve.append(- ((b - a) / 2 ) * A_I_J[i][j] * K[i][j] / koef)
        SOLVE.append(solve)
    SOLVE = np.array(SOLVE)


    # уравнение типа SOLVE * Y = F готово
    # Y = SOLVE **(-1) * F
    Y = np.dot(np.linalg.inv(SOLVE), f_values)

    return nodes_array, Y
if __name__ == '__main__':
    n = 6
    F = I_FUNC(-1, 1)
    x_array, y_array = solve(F, n)
    pol = LAGR.create_Lagrange_polynomial(x_array, y_array)

    x = np.arange(-2, 2, 0.0001)


    plt.plot(x_array, y_array)
    plt.title("n = " + str(n))
    plt.grid()
    plt.scatter(x_array, y_array)
    plt.show()



