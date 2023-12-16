def Simpson_method(range_1, range_2, func, delta_x = 0.001):
    n = int((range_2 - range_1)/delta_x) # кол-во разбиений
    sum = 0

    for i in range(n):
        try:
            func_result = func(range_1 + i*delta_x) + 4 * func(range_1 + (i + 0.5)*delta_x)
            func_result += func(range_1 + (i + 1)*delta_x)
            sum += func_result*delta_x/6
        except:
            sum += 0
    return sum
def trapezoid_method(range_1, range_2, func, delta_x = 0.001):
    n = int((range_2 - range_1)/delta_x) # кол-во разбиений
    sum = 0

    for i in range(n):
        try:
            sum += ((func(range_1 + i*delta_x) + func(range_1 + (i+1)*delta_x))/2)*delta_x
        except:
            sum += 0
    return sum
