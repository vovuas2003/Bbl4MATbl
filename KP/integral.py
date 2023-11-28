def trap(f, a, b, h):
    import numpy as np
    x = np.arange(a, b + h, h)
    su = 0
    for i in range(len(x) - 1):
        su += f(x[i]) + f(x[i + 1])
    return su * h / 2

def trap_auto(f, a, b, eps):
    h = 1
    i1 = trap(f, a, b, h)
    i2 = trap(f, a, b, h / 2)
    while(abs(i1 - i2) > (eps * (2**2 - 1))):
        h /= 2
        i1 = i2
        i2 = trap(f, a, b, h / 2)
    print("I = " + str(i2))
    print("h = " + str(h / 2))

def trap_auto_start(f, a, b, eps, h):
    i1 = trap(f, a, b, h)
    i2 = trap(f, a, b, h / 2)
    while(abs(i1 - i2) > (eps * (2**2 - 1))):
        h /= 2
        i1 = i2
        i2 = trap(f, a, b, h / 2)
    print("I = " + str(i1))
    print("h = " + str(h))

def simp(f, a, b, h):
    import numpy as np
    x = np.arange(a, b + h, h)
    su = 0
    for i in range(len(x) - 1):
        su += f(x[i]) + f(x[i + 1]) + 4 * f(x[i] + h / 2)
    return su * h / 6

def simp_auto(f, a, b, eps):
    h = 1
    i1 = simp(f, a, b, h)
    i2 = simp(f, a, b, h / 2)
    while(abs(i1 - i2) > (eps * (2**4 - 1))):
        h /= 2
        i1 = i2
        i2 = simp(f, a, b, h / 2)
    print("I = " + str(i2))
    print("h = " + str(h / 2))

def simp_auto_start(f, a, b, eps, h):
    i1 = simp(f, a, b, h)
    i2 = simp(f, a, b, h / 2)
    while(abs(i1 - i2) > (eps * (2**4 - 1))):
        h /= 2
        i1 = i2
        i2 = simp(f, a, b, h / 2)
    print("I = " + str(i1))
    print("h = " + str(h))

def gauss2(f, a, b):
    t = 1 / (3**0.5)
    p = (b - a) / 2
    f1 = f((b + a) / 2 + p * t)
    f2 = f((b + a) / 2 - p * t)
    return p * (f1 + f2)

def gauss3(f, a, b):
    t = (3 / 5)**0.5
    p = (b - a) / 2
    f1 = f((b + a) / 2 + p * t)
    f2 = f((b + a) / 2)
    f3 = f((b + a) / 2 - p * t)
    return p * (5 * (f1 + f3) / 9 + 8 * f2 / 9)
