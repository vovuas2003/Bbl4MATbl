def getF(k : int, n : int, ArgValues : list, FunctionValues : list) -> float:
    # Находит разделенную разность F(tk, ... tn)
    # ArgValues      - массив из n         значений аргумента [t0, ... tn]
    # FunctionValues - массив из n - k + 1 значений функции в точках [tk, ... tn]
    if (k == n):
        return FunctionValues[0]
    F2 = getF(k + 1, n, ArgValues, FunctionValues[1:])
    F1 = getF(k, n - 1, ArgValues, FunctionValues[:-1])
    t2 = ArgValues[n]
    t1 = ArgValues[k]
    DivDiff = (F2 - F1) / (t2 - t1)
    return DivDiff

def getNewtonPolinom(ArgValues : list, FunctionValues : list) -> list:
    # Возвращает коэффициенты интерполяционного полинома Ньютона (разделенные разности)
    n = len(ArgValues)
    DivDiffs = []
    for k in range(n):
        F = getF(0, k, ArgValues, FunctionValues)
        DivDiffs.append(F)
    return DivDiffs
    
def getNewtonValue(t : float, ArgValues : list, NewtonPol : list) -> float:
    # Возвращает значение полинома Ньютона в точке t
    Result = 0
    n = len(ArgValues)
    for k in range(n):
        Mult = 1
        for i in range(k):
            Mult *= (t - ArgValues[i])
        Result += Mult * NewtonPol[k]
    return Result

def getNewtonValues(Args : list, ArgValues : list, FunctionValues : list) -> list:
    # Возвращает масcив значений полинома Ньютона для каждого значения аргумента из Args
    # ArgValues      - массив из n значений аргумента t -- [t0, ... tn]
    # FunctionValues - массив из n значений функции в точках [t0, ... tn]
    NewtonValues = []
    NewtonPol = getNewtonPolinom(ArgValues, FunctionValues)
    for Arg in Args:
        NewtonValues.append(getNewtonValue(Arg, ArgValues, NewtonPol))
    return NewtonValues

def getChebZeros(Start : float, Stop : float, N : int) -> float:
    # Возвращает массив из N нулей полинома Чебышева на отрезке [Start, Stop]
    Zeros = []
    HalfSum = (Start + Stop) / 2.0
    HalfDiff = (Stop - Start) / 2.0
    for ZeroNum in range(1, N + 1):
        Zero = HalfSum + HalfDiff * math.cos((2 * ZeroNum - 1) * math.pi / (2 * N))
        Zeros.append(Zero)
    return Zeros

def main():
    StartSegment = -1
    EndSegment = 1
    
    plt.figure(figsize = (10, 10))
    plt.title("Интерполяционный многочлен Ньютона с узлами в нулях полинома Чебышева при разных n")
    
    Args = np.arange(-1, 1.01, 0.01)
    Vals = [Function1(Arg) for Arg in Args]
    plt.plot(Args, Vals, 'r', label = "График исходной функции")
    
    ArrayN = [4, 6, 10]
    
    for N in ArrayN:
        NewtonArgValues = getChebZeros(StartSegment, EndSegment, N)
        FunctionValues = []
        for Arg in NewtonArgValues:
            FunctionValues.append(Function1(Arg))
        plt.scatter(NewtonArgValues, FunctionValues, marker = "^")
        Args = np.arange(-1, 1.01, 0.01)
        NewtonVals = getNewtonValues(Args, NewtonArgValues, FunctionValues)
        NameGraph = "Ньютон n =" + str(N)
        plt.plot(Args, NewtonVals, label = NameGraph)

    plt.grid()
    plt.legend()
    plt.savefig("NewtonWithCheb.png")
    plt.show()


if __name__ == '__main__':
    main()
