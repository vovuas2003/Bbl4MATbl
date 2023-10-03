import matplotlib.pyplot as plt
import numpy as np

N = 15

#1-я норма вектора
def norm_v(x):
    return max(abs(x))

def make_A(n, a):
    A = np.eye(n)
    A *= 2
    for i in range(n - 1):
        A[i, i + 1] = -1 - a
    for i in range(1, n):
        A[i, i - 1] = -1 + a
    return A
def make_f(n, a):
    f = np.zeros(n)
    f[0] = 1 - a
    f[-1] = 1 + a
    return f
def LDU(A):
    n = len(A)
    L = np.zeros((n, n))
    D = np.zeros((n, n))
    U = np.zeros((n, n))
    for i in range(n):
        D[i][i] = A[i][i]
    for i in range(1, n):
        for j in range(i):
            L[i][j] = A[i][j]
    for i in range(n - 1):
        for j in range(i + 1, n):
            U[i][j] = A[i][j]
    return L, D, U

def solve(A, f):
    n = 1
    L, D, U = LDU(A)
    prev = np.asarray([f[i]/A[i][i] for i in range(len(f))])
    #prev = np.zeros(len(f))
    x = -1 * (np.linalg.inv(L + D).dot(U)).dot(prev) + (np.linalg.inv(L + D)).dot(f)
    while(norm_v(x - prev) > 10**-6):
        n +=1
        prev = x
        x = -1 * (np.linalg.inv(L + D).dot(U)).dot(prev) + (np.linalg.inv(L + D)).dot(f)
    return x, n

x=[]
y=[]
a = [0.01*i for i in range(101)]
for i in range(len(a)):
    _, n = solve(make_A(N, a[i]), make_f(N, a[i]))
    x.append(a[i])
    y.append(n)


plt.plot(x, y)
plt.title("n = " + str(N))
plt.xlabel("alpha")
plt.ylabel("число итераций")
plt.tight_layout()
plt.show()
