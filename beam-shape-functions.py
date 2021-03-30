from sympy import *

x = S("x")
A = S("A")
B = S("B")
C = S("C")
D = S("D")

L = S("L")

v_1 = S("v_1")
v_2 = S("v_2")
v_4 = S("v_4")
v_5 = S("v_5")

uy = A*x**3 + B*x**2 + C*x + D
duy_dx = uy.diff(x)


#condiciones de borde:

#    uy(0) - v_1 = 0
cb1 = uy.subs(x, 0) - v_1
cb2 = duy_dx.subs(x, 0) - v_2
cb3 = uy.subs(x, L) - v_4
cb4 = duy_dx.subs(x, L) - v_5


print(cb1)
print(cb2)
print(cb3)
print(cb4)

sol = solve((cb1, cb2, cb3, cb4), (A, B, C, D))

A_sol = sol[A]
B_sol = sol[B]
C_sol = sol[C]
D_sol = sol[D]

print("A = ", A_sol)
print("B = ", B_sol)
print("C = ", C_sol)
print("D = ", D_sol)

uy_sol = uy.subs(A, A_sol).subs(B, B_sol).subs(C, C_sol).subs(D, D_sol)

print(uy_sol)



N1 = uy_sol.subs(v_1,1).subs(v_2,0).subs(v_4,0).subs(v_5,0)
N2 = uy_sol.subs(v_1,0).subs(v_2,1).subs(v_4,0).subs(v_5,0)
N4 = uy_sol.subs(v_1,0).subs(v_2,0).subs(v_4,1).subs(v_5,0)
N5 = uy_sol.subs(v_1,0).subs(v_2,0).subs(v_4,0).subs(v_5,1)

print("N1(x) = ", N1)
print("N2(x) = ", N2)
print("N4(x) = ", N4)
print("N5(x) = ", N5)