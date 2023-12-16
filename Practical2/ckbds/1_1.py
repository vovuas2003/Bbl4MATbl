import numpy as np
import math
import matplotlib.pyplot as plt

def Function1(X : float) -> float:
    # Функция из задачи 1
    return 1.0 / (1 + 25 * X * X)

def getl(t : float, k : int, ArgValues : list) -> float:
    # Находит базисную функцию Лагранжа lk(t)
    # ArgValues - массив из n значений аргумента [t0, ... tn]
    n = len(ArgValues)
    lk = 1
    for j in range(n):
        Denom = ArgValues[k] - ArgValues[j]
        lk *= ((t - ArgValues[j]) / Denom) if k != j else 1
    return lk

def getLagrangePolinom(t : float, ArgValues : list, FunctionValues : list) -> float:
    # Возвращает значение полинома Лагранжа в точке t
    n = len(ArgValues)
    Value = 0
    for k in range(n):
        Value += (getl(t, k, ArgValues) * FunctionValues[k])
    return Value

def getLagrangeValues(Args : list, ArgValues : list, FunctionValues : list) -> list:
    # Возвращает массив значений полинома Лагранжа для каждого значения аргумента из Args
    # ArgValues      - массив из n значений аргумента t -- [t0, ... tn]
    # FunctionValues - массив из n значений функции в точках [t0, ... tn]
    Values = []
    for Arg in Args:
        Values.append(getLagrangePolinom(Arg, ArgValues, FunctionValues))
    return Values

def divideSegment(Start : float, Stop : float, N : int) -> float:
    # Возвращает массив из N равноотстоящих точек отрезка [Start, Stop]
    Points = []
    Step = (Stop - Start) / (N - 1)
    for PointNum in range(N):
        Points.append(Start + PointNum * Step)
    return Points

def main():
    StartSegment = -1
    EndSegment = 1
    
    plt.figure(figsize = (10, 10))
    plt.title("Интерполяционный многочлен Лагранжа при разных n")
    
    Args = np.arange(-1, 1.01, 0.01)
    Vals = [Function1(Arg) for Arg in Args]
    plt.plot(Args, Vals, 'r', label = "График исходной функции")
    
    ArrayN = [4, 6, 10]
    
    for N in ArrayN:
        ArgValues = divideSegment(StartSegment, EndSegment, N)
        FunctionValues = []
        for Arg in ArgValues:
            FunctionValues.append(Function1(Arg))
        plt.scatter(ArgValues, FunctionValues, marker = "^")
        Args = np.arange(-1, 1.01, 0.01)
        LagrangeVals = getLagrangeValues(Args, ArgValues, FunctionValues)
        NameGraph = "Лагранж n =" + str(N)
        plt.plot(Args, LagrangeVals, label = NameGraph)

    plt.grid()
    plt.legend()
    plt.savefig("Lagrang.png")
    plt.show()

if __name__ == '__main__':
    main()
