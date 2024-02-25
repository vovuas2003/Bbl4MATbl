import numpy as np
import matplotlib.pyplot as plt

def Jacobian( f, x ):

  h = 1.0e-4

  J = np.zeros( ( x.size, x.size ) )

  fn = f(x)

  for i in np.arange( 0, x.size, 1 ):

    x_old = np.copy( x )

    x_old[ i ] = x_old[ i ] + h

    fn_1 = f(x_old)

    J[:,i] = ( fn_1 - fn ) / h

  return J, fn

def NewtonMethod( f, x, eps = 1.0e-6 ):

  max_iter = 10000

  for i in np.arange( 0, max_iter, 1 ):

    J, fn = Jacobian( f, x )

    if np.sqrt( np.dot( fn, fn ) / x.size ) < eps:

      return x, i

    dx = np.linalg.solve( J, fn )

    x = x - dx

def Iterarion( f, y0, tBEG, tEND, tau, alpha ):

  def F(y_next):

    return y_next - tau * alpha * f(t[i], y_next) - y[i] - tau * ( 1. - alpha) * f( t[i], y[i] )

  t = np.arange( tBEG, tEND, tau )

  y = np.zeros( ( t.size, 2 ) )

  y[0] = y0

  for i in np.arange(0, t.size-1, 1):

    y_next = y[ i ] + tau * f( t[ i ], y[ i ] )

    y[ i + 1 ], iter = NewtonMethod( F, y_next )

  return t, y

def f(t, y):

  f_fun = np.zeros( 2 )

  f_fun[ 0 ] = y[ 0 ] - y[ 0 ] * y[ 1 ]

  f_fun[ 1 ] = -y[ 1 ] + y[ 0 ] * y[ 1 ]

  return f_fun

tBEG = 0.

tEND = 10.

tau = 0.001

y0 = np.asarray([ 2., 2. ])

alpha = 0.5

t, y = Iterarion( f, y0, tBEG, tEND, tau, alpha )

for n in np.arange(0,2,1):

  r = y[:,n]

  label = r'$y_1$'

  mark = '-'

  if n == 1:

    label = r'$y_2$'

    mark = '--'

  plt.plot( t, r, mark, label = label )

plt.legend()
plt.xlabel(r'$t$')
plt.grid()
plt.show()
