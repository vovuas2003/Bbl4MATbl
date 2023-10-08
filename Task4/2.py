eps = 4 #10^-eps

def f1(x, y):
    return x * y + x**2 - 1.03
def f2(x, y):
    return -2 * x**3 + y**2 - 1.98

xp = 1
yp = 2
x = 1
y = 2
e = 10**-15
def my_1_3(x):
    n = 1
    if(x < 0):
        x *= -1
        n *= -1
    return n * (x ** (1/3))
while abs(f1(x, y)) > e or abs(f2(x, y)) > e:
    xp = x
    yp = y
    x = my_1_3((yp**2 - 1.98) / 2)
    y = -1 * xp + 1.03 / xp
print("x = " + str(x) +  ", y = " + str(y))

xr = x
yr = y
xp = 1
yp = 2
x = 1
y = 2
e = 10**-eps
n = 0
while abs(x - xr) > e or abs(y - yr) > e:
    n += 1
    xp = x
    yp = y
    x = my_1_3((yp**2 - 1.98) / 2)
    y = -1 * xp + 1.03 / xp
print(str(n) + " steps for eps = 10^-" + str(eps))
