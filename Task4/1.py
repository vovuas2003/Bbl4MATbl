import matplotlib.pyplot as plt
import numpy as np
import math

st = 6 #погрешность 10^-st

def f1(x):
    return math.sqrt(1 - x**2)

def f2(x):
    return -math.sqrt(1 - x**2)

def f3(x):
    return math.tan(x)

x = np.linspace(-1, 1, 2 * 10**st + 1)
y1 = [f1(i) for i in x]
y2 = [f2(i) for i in x]
y3 = [f3(i) for i in x]
rx = []
e = 10**-st
for i in x:
    if abs(f1(i) - f3(i)) < e or abs(f2(i) - f3(i)) < e:
        rx.append(i)
ry = [f3(i) for i in rx]

plt.figure(figsize = (8, 8))
plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)
plt.scatter(rx, ry, color = 'r', s = 100)
for a, b in zip(rx, ry): 
    plt.text(a + 0.1, b, str(round(a, 2)) + "; " + str(round(b, 2)), size = 'large')
plt.plot(x, y1, 'b', zorder = -1)
plt.plot(x, y2, 'b', zorder = -1)
plt.plot(x, y3, 'g', zorder = -1)
#plt.grid()
plt.tight_layout()
plt.show()
