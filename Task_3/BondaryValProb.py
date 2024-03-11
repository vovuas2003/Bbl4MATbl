import numpy as np
import matplotlib.pyplot as plt

def Lam( k ):
    global h, X
    return -4. / h ** 2 * np.sin( np.pi * k * h / 2. / X ) ** 2

def Omega( k, x ):
    global X
    return np.sqrt( 2. / X ) * np.sin( np.pi * k * x / X )

def g( x ):
    return x ** 3

def MethTr(f, x):
    global h
    S = 0.
    for i in np.arange(0, x.size - 1, 1 ):
        S += h / 2. * ( f[ i ] + f[ i + 1 ] )
    return S

def Matrix( x ):
    global h, X
    global Omega, MethTr
    A = np.zeros( ( x.size-1, x.size-1 ) )
    for i in np.arange( 0, x.size-1, 1 ):
        for j in np.arange( 0, x.size-1, 1 ):
            A[ i ][ j ] = MethTr( Omega( i + 1, x ) * Omega( j + 1, x ), x )
    return A

def F( x ):
    global g, Omega
    f = np.zeros( x.size-1 )
    for i in np.arange( 0, x.size-1, 1 ):
        f[ i ] = MethTr( g( x ) * Omega( i + 1, x ) )
    return f

def Cg( A, f ):
    return np.linalg.solve( np.linalg.inv( A ), f )

h = 0.001
X = 1.0
x1 = np.arange( 0, X+h, h )
G = np.zeros( x1.size - 1 )
A = Matrix( x1 )
Fg = F( x1 )
C = Cg(A, Fg)

for i in np.arange( 0, x1.size-1, 1 ):
    for j in np.arange( 0, x1.size-1, 1 ):
        G[ i ] += C[ j ] * Omega( j + 1, x1[ i ] )

plt.figure( figsize = ( 10, 6 ) )
plt.title("Close plot to continue")
plt.rc('font', **{'size' : 15})
plt.plot( x1, g(x1), label = r'$g(x)$')
plt.plot(x1[:-1],G, label = r'$G(x)$')
plt.ylabel(r'$g(x)$')
plt.xlabel(r'$x$')
plt.legend()
plt.grid()
plt.show()

def C_coeff( Cg ):
    global Lam
    C = np.zeros( Cg.size-1 )
    for i in np.arange( 0, Cg.size-1, 1 ):
        C[ i ] = Cg[ i ]  / Lam( i + 1 )
    return C

def U( A, f, x ):
    global C_coeff, Cg, Omega
    C = C_coeff( Cg( A, f ) )
    u = np.zeros( x.size )
    for i in np.arange( 0, x.size-1, 1 ):
        for j in np.arange( 0, x.size-1, 1 ):
            u[ i ] += C[ j ] * Omega( j + 1, x[ i ] )
    return u

u1 = U( A, Fg, x1 )
plt.figure( figsize = ( 10, 6 ) )
plt.title("Close plot to continue")
plt.rc('font', **{'size' : 15})
plt.plot( x1, u1 )
plt.ylabel(r'$u(x)$')
plt.xlabel(r'$x$')
plt.grid()
plt.show()

def f(x, y, v):
    global g
    return g( x )

def MethEl( x, y, v ):
    global h
    v = f( x, y, v ) * h + v
    y = v * h + y
    return v, y

def Alpha( p, alpha0, y, phi ):
    r = y - phi
    return alpha0 - 1.0 / p * r

def CalYstrel( x, y, Ylast, Scheme):
    global Alpha
    N = x.shape[ 0 ]
    v = np.zeros( N )
    v[ 0 ] = 1.13
    y_st = np.copy( y )
    v_st = np.copy( v )
    v_st[0] = v[0] + 0.005
    while np.fabs( y[ -1 ] - Ylast ) > 1.0e-2:
        for i in np.arange( 0, N - 1, 1 ):
            v[ i + 1 ], y[ i + 1 ] = Scheme( x[ i ], y[ i ], v[ i ] )
            v_st[ i + 1 ], y_st[ i + 1 ] = Scheme( x[ i ], y_st[ i ], v_st[ i ] )
        v[ 0 ] = Alpha( ( y_st[-1] - y[-1] ) / 0.005, v[ 0 ], y[ -1 ], Ylast )
        v_st[0] = v[0] + 0.005
    print( 'Last alpha =', v[ 0 ] )
    return y

h = 1.0e-3
x = np.arange( 0, 1, h )
y = np.copy( x )
y[ 0 ] = 0.
y = CalYstrel( x, y, 0, MethEl )
plt.figure( figsize = ( 10, 6 ) )
plt.title("Close plot to continue")
plt.rc('font', **{'size' : 15})
plt.plot( x, y )
plt.ylabel(r'$u(x)$')
plt.xlabel(r'$x$')
plt.grid()
plt.show()

plt.figure( figsize = ( 10, 6 ) )
plt.title("Last plot")
plt.rc('font', **{'size' : 15})
plt.plot( x, y, label = 'Метод стрельбы' )
plt.plot( x1, u1, label = 'Метод Фурье' )
plt.ylabel(r'$u(x)$')
plt.xlabel(r'$x$')
plt.legend()
plt.grid()
plt.show()
