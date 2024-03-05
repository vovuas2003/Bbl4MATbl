import numpy as np
import matplotlib.pyplot as plt

N = 3
def f(t, y):
    f_fun = np.zeros(N)
    f_fun[0] = y[0] * (1 - 0.5 * y[0] - (2 * y[1]) / (7 * ((y[2])**2)))
    f_fun[1] = y[0] * (2 * y[2] - 3.5 * y[0] * ((y[2])**2) - 0.5 * y[1])
    f_fun[2] = (2 - 7 * y[2] * y[0]) / 100
    return f_fun
LABELS = [r'$x$', r'$y$', r'$a$'] # r'$y_1$', r'$y_2$'
MARKS = [] # '-', '--'
XLABEL = r'$t$'

def main():
    tBEG = 0.
    tEND = 0.5
    tau = 0.001
    y0 = np.asarray([1.5, 10., 0.0005])
    if len(y0) != N:
        print("Len init values error!")
        return
    alpha = 0.5
    t, y = Iterarion(f, y0, tBEG, tEND, tau, alpha)
    plt.figure(figsize = (13.5, 6.3))
    if len(LABELS) == 0 or len(LABELS) != N:
            print("Use default labels.")
    for n in np.arange(0, N, 1):
        r = y[:, n]
        if len(LABELS) == 0 or len(LABELS) != N:
            label = r'$y_' + str(n + 1) + '$'
        else:
            label = LABELS[n]
        if len(MARKS) == 0:
            mark = '-' * ((n % 2) + 1)
        else:
            mark = MARKS[n % len(MARKS)]
        plt.plot(t, r, mark, label = label)
    plt.legend()
    plt.xlabel(XLABEL)
    plt.tight_layout()
    plt.grid()
    plt.show()

def Jacobian(f, x):
    h = 1.0e-4
    J = np.zeros((x.size, x.size))
    fn = f(x)
    for i in np.arange(0, x.size, 1):
        x_old = np.copy(x)
        x_old[i] = x_old[i] + h
        fn_1 = f(x_old)
        J[:, i] = (fn_1 - fn) / h
    return J, fn

def NewtonMethod(f, x, eps = 1.0e-6):
    max_iter = 10000
    for i in np.arange(0, max_iter, 1):
        J, fn = Jacobian(f, x)
        if np.sqrt(np.dot(fn, fn) / x.size) < eps:
            return x, i
        dx = np.linalg.solve(J, fn)
        x = x - dx

def Iterarion(f, y0, tBEG, tEND, tau, alpha):
    def F(y_next):
        return y_next - tau * alpha * f(t[i], y_next) - y[i] - tau * (1. - alpha) * f(t[i], y[i])
        #return y_next - tau * alpha * f(t[i + 1], y_next) - y[i] - tau * (1. - alpha) * f(t[i], y[i])
    t = np.arange(tBEG, tEND, tau)
    y = np.zeros((t.size, N))
    y[0] = y0
    for i in np.arange(0, t.size-1, 1):
        y_next = y[i] + tau * f(t[i], y[i])
        y[i + 1], iter = NewtonMethod(F, y_next)
    return t, y

if __name__ == "__main__":
    main()
