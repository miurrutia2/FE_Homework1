from sympy import *
from sympy.matrices import Matrix

x = S("x")
L = S("L")

N1 = 1- 3*x**2/L**2 + 2*x**3/L**3
N2 = x -2*x**2/L + x**3/L**2
N3 = 3*x**2/L**2 - 2*x**3/L**3
N4 = -x**2/L + x**3/L**2

B1 = N1.diff(x,2)
B2 = N2.diff(x,2)
B3 = N3.diff(x,2)
B4 = N4.diff(x,2)

print(B1)
print(B2)
print(B3)
print(B4)

B = Matrix([[B1, B2, B3, B4]])

#print(B)
integral  = integrate(B.T*B, (x, 0, L))

#print(pretty(integral))
