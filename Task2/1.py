#Ответ: при n = 3.
import matplotlib.pyplot as plt
import math

def u(t):
    return math.sin(t)
t = 0.5

#факториал
def fac(x):
    p = 1
    while x != 0:
        p *= x
        x -=1
    return p
#число сочетаний из n по k
def soch(n, k):
    return fac(n) // (fac(n - k) * fac(k))
#n-я проивзодная численно
def pr(f, x, h, n):
    s = 0
    for k in range(0, n + 1):
        s += soch(n, k) * f(x + k * h) * ((-1) ** (n - k))
    return s / (h ** n)
#n-я частичная сумма ряда Маклорена
def macl(f, x, h, n):
    s = 0
    for i in range(0, n + 1):
        s += (pr(f, 0, h, i) * (x ** i)) / fac(i)
    return s
#настройки
h = 10**-3
n = 6

x = [i for i in range(1, n + 1)]
y = [macl(u, t, h, i) for i in x]
ut = u(t)
y = [abs(ut - u) for u in y]
x1 = [1, n]
y1 = [10**-3, 10**-3]
plt.scatter(x, y)
plt.plot(x1, y1, 'r')
plt.show()
