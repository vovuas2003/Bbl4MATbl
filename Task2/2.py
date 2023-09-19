import matplotlib.pyplot as plt
import numpy as np
import math
n = 5
a = np.eye(n)
a *= -2
for i in range(n - 1):
    a[i, i + 1] = 1
for i in range(1, n):
    a[i, i - 1] = 1
a_1 = np.linalg.inv(a) #обратная матрица
def n1(a):
    return len(a)
print(n1(a))
aa = np.dot(np.conj(a), a)
l, _ = np.linalg.eig(aa) #массив собственных чисел
l = [abs(i) for i in l]
print(max(l))
n3 = math.sqrt(max(l)) #3-я норма матрицы
