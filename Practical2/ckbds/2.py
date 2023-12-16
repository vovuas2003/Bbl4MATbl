import scipy
from scipy import integrate

def Function2(X : float) -> float:
    # Функция из задания 2
    return math.log(100.0 - X) / (10.0 - math.sqrt(X))

def calculateSimpsonByPoints(Functions : list, a : float, b : float) -> float:
    # Считает интеграл от поточечно заданной функции со значениями Functions[i] на отрезке [a, b] методом Симпсона
    Int = 0
    N = len(Functions)
    k = int(N / 2)
    h = (b - a) / N
    for i in range(1, k):
        F1 = Functions[2*i]
        F2 = Functions[2*i - 1]
        F3 = Functions[2*i - 2]
        Int += h / 3.0 * (F1 + 4 * F2 + F3)
    return Int

def getNextNewtonIteration(Xk : float, P : float, P1 : float) -> float:
    # Возвращает следующий член в итерационном методе Ньютона
    # Процесс представлен следующим итерационным соотношением: X_{k+1} = X_{k} - P(Xk) / Р1(Xk)
    # P - значение полинома Лежандра для Xk
    # P1 - значение произвожной полинома Лежандра для Xk
    return Xk - P / P1

def getLejanPol(X : float, N : int) -> float:
    # Возвращает N-ый полином Лежандра для X
    if (N == 0):
        return 1
    if (N == 1):
        return X
    return (2.0 * N + 1) * X * getLejanPol(X, N - 1) / (N + 1) - N * getLejanPol(X, N - 2) / (N + 1)

def getLejanDerr(X : float, N : int) -> float:
    # Возвращает первую производную N-ного полинома Лежандра для X
    return N * (getLejanPol(X, N - 1) - X * getLejanPol(X, N)) / (1 - X * X)

def getLejanZeros(N : int) -> list:
    # Возвращает N нулей полинома Лежандра
    #, вычисленные итеративно по методу Ньютона с начальным приближением: X0 = cos(pi(4i - 1)/(4N + 2))
    Zeros = []
    Epsilon = 1e-3
    for i in range(1, N + 1):
        Xk = math.cos(math.pi * (4 * i - 1) / (4 * N + 2))
        Xk1 = getNextNewtonIteration(Xk, getLejanPol(Xk, N), getLejanDerr(Xk, N))
        while (abs(Xk - Xk1) > Epsilon):
            Xk = Xk1
            Xk1 = getNextNewtonIteration(Xk, getLejanPol(Xk, N), getLejanDerr(Xk, N))
        Zeros.append(Xk1)
    return Zeros

def changeVars(Start : float, Stop : float, Vars : list) -> list:
    # Делает замену переменных в узлах квадратуры для перехода от интегрирования
    # по [-1, 1] к интегрированию по [a, b]
    HalfSum = (Start + Stop) / 2.0
    HalfDiff = (Stop - Start) / 2.0
    NewVars = [HalfSum + HalfDiff * T for T in Vars]
    return NewVars

def calculateWithGaussQuadrature(F : "function", a : float, b : float, N : int) -> float:
    # Считает интеграл от F на отрезке [a, b] методом квадратур Гаусса по N узлам
    Int = 0
    NodesT = getLejanZeros(N)
    NodesX = changeVars(a + 1e-1, b, NodesT)
    for k in range(1, N + 1):
        # Интегрирование методом Симпсона по 1000 точкам
        Args = np.arange(a, b, 0.001)
        BaseLagranValues = [getl(i, k - 1, NodesX) for i in Args]
        Ck = calculateSimpsonByPoints(BaseLagranValues, a, b)
        Fk = F(NodesX[k - 1])
        Int += Ck * Fk
    return Int
    

def main():
    StartSegment = 0
    EndSegment = 10
    
    ArrayN = np.arange(2, 10)
    Errors =[]
    
    Exact = scipy.integrate.quad(Function2, StartSegment, EndSegment)
    print("Точное решение и его ошибка (I, delta I) =", Exact)
    
    for N in ArrayN:
        Gauss = calculateWithGaussQuadrature(Function2, StartSegment, EndSegment, N)
        print("\nМетод квадратур Гаусса при N =", N, ": I =", format(Gauss, '.10f'))
        Errors.append(abs(Gauss - Exact[0]) * 100 / Exact[0])

    plt.figure(figsize = (10, 10))
    plt.title("Кривая зависимости относительной ошибки интегрирования от количества узлов,\n интерполированная полиномом Лагранжа:)")
    plt.scatter(ArrayN, Errors, marker = '^', color = 'r')
    Args = np.arange(min(ArrayN), max(ArrayN), 0.01)
    NewtonVals = getLagrangeValues(Args, ArrayN, Errors)
    plt.plot(Args, NewtonVals, color = 'g')
    plt.xlabel('Количество узлов n')
    plt.ylabel('Ошибка, %')
    plt.grid()
    plt.savefig("Gauss.png")
    plt.show()


if __name__ == '__main__':
    main()
