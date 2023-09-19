import matplotlib.pyplot as plt
import numpy as np
import math

N = 100

#1-я норма матрицы
def n1(a):
    a = abs(a)
    su = a.sum(axis = 1) #массив сумм строк
    return max(su)
#2-я норма матрицы
def n2(a):
    a = abs(a)
    su = a.sum(axis = 0) #массив сумм столбцов
    return max(su)
#3-я норма матрицы
def n3(a):
    aa = np.dot(np.transpose(np.conj(a)), a)
    l, _ = np.linalg.eig(aa) #массив собственных чисел
    l = [abs(i) for i in l]
    return math.sqrt(max(l))

x = []
y1 = []
y2 = []
y3 = []
for n in range(1, N + 1):
    a = np.eye(n)
    a *= -2
    for i in range(n - 1):
        a[i, i + 1] = 1
    for i in range(1, n):
        a[i, i - 1] = 1
    a_1 = np.linalg.inv(a) #обратная матрица
    x.append(n)
    y1.append(n1(a) * n1(a_1))
    y2.append(n2(a) * n2(a_1))
    y3.append(n3(a) * n3(a_1))

plt.figure()
plt.subplot(1, 3, 1)
plt.plot(x, y1)
plt.tight_layout()
plt.subplot(1, 3, 2)
plt.plot(x, y2)
plt.tight_layout()
plt.subplot(1, 3, 3)
plt.plot(x, y3)
plt.tight_layout()
plt.show()
