import math
import matplotlib.pyplot as plt
import numpy as np
import Lagrange_method as Inter
import Integral as INTEGRAL
import scipy.integrate as sp

def main_func(x):
    return  math.log(100 - x)/(10 - math.sqrt(x))

def grafic(a, b, func, x_array = [], y_array = [], delta_x = 0.001):

    x = list(np.arange(a, b, delta_x))
    y = []
    for i in x:
        y.append(func(i))

    plt.plot(x, y)
    plt.grid()
    plt.scatter(x_array, y_array)

def calculate_the_integral(a, b, func, n):

    # работаем с полиномом лежандра
    def polinom_Legandre(x, n):
        if n == 0:
            return 1
        elif n == 1:
            return x
        elif n > 1:
            return ((2*n - 1)/n) * x * polinom_Legandre(x, n - 1) - ((n - 1)/n)*polinom_Legandre(x, n - 2)
        else:
            print("Ошибка: такого n быть не может")
    def derivative_polinom_Legandre(x, n):
        return (n/(1 - x**2)) * (polinom_Legandre(x, n - 1) - x*polinom_Legandre(x, n))
    def Zero_of_Legendre_polynomials(k, n, eps = 0.0001):
        x_0 = math.cos( (math.pi*(4*k - 1)) / (4*n + 2) )
        x_k = x_0 - polinom_Legandre(x_0, n)/derivative_polinom_Legandre(x_0, n)

        while abs(x_k - x_0) > eps:
            x_0 = x_k
            x_k = x_0 - polinom_Legandre(x_0, n)/derivative_polinom_Legandre(x_0, n)


        return x_k


    # задаем узлы и их значения
    x_array = []

    for k in range(1, n + 1):
        x_array.append(Zero_of_Legendre_polynomials(k, n))



    # создаем базовые многочлены многочленов Лагранжа

    basic_polinom = []
    for i in range(len(x_array)):
        basic_polinom.append(Inter.create_basic_polynomial(x_array, i))



    # создаем веса квадратурной формулы

    weight_formula = []
    for i in range(len(basic_polinom)):
        weight_formula.append(INTEGRAL.Simpson_method(-1, 1, basic_polinom[i]))

    # вычисляем интеграл
    sum = 0
    for i in range(len(weight_formula)):
        sum += func((b+a)/2 + x_array[i]*(b-a)/2 )*weight_formula[i]
    sum = sum*(b-a)/2


    return sum




if __name__ == '__main__':

    #grafic(0, 10, main_func)
    #plt.show()
    a = 0
    b = 10

    I_exact, I_mistake = sp.quad(main_func, a, b)
    print("Точное значение интеграла: " + str(I_exact))

    #строим график зависимости ошибки от кол-ва узлов
    x_array = [i for i in range(2, 15)]
    y_array = [abs(calculate_the_integral(a, b, main_func, n) - I_exact) for n in x_array]
    plt.plot(x_array, y_array)
    plt.grid()
    plt.show()


