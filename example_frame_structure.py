from numpy import array, pi, zeros, ix_
from beam_element import beam_element
import matplotlib.pylab as plt
from scipy.linalg import solve


xy = array([
	[0,0],
	[4,3],
	[9,3],
	[9,0]
	])

conec = array([
	[0,1],
	[1,2],
	[2,3]
	], dtype = int)

t = 20e-3
r = 400e-3

properties0 = {}
properties1 = {}

properties0["E"] = 200e9
properties0["A"] = pi*(r**2 - (r-t)**2)
properties0["I"] = pi*(r**4/4 - (r-t)**4/4)

properties0["qx"] = 0
properties0["qy"] = 0

properties1["E"] = 200e9
properties1["A"] = pi*(r**2 - (r-t)**2)
properties1["I"] = pi*(r**4/4 - (r-t)**4/4)

properties1["qx"] = 0
properties1["qy"] = -5000.

properties = [properties0, properties1, properties0]

Nnodes = xy.shape[0]
Nelems = conec.shape[0]

NDOFs_per_node = 3
NDOFs = 3*Nnodes

K = zeros((NDOFs, NDOFs))
f = zeros((NDOFs, 1))

for e in range(Nelems):
	ni = conec[e,0]
	nj = conec[e,1]

	#print(f"e = {e} ni = {ni} nj = {nj}")

	xy_e = xy[[ni, nj], :]
	#print(f"xy_e = {xy_e}")

	ke, fe = beam_element(xy_e, properties[e])


	#print(f"ke = {ke}")
	d = [3*ni, 3*ni+1, 3*ni+2, 3*nj, 3*nj+1, 3*nj+2]

	for i in range(2*NDOFs_per_node):
		p = d[i]
		for j in range(2*NDOFs_per_node):
			q = d[j]
			K[p,q] += ke[i,j]
		f[p] += fe[i]

q =5000.
L = 5.0

f[4] = -q*L/2
f[5] = -q*L**2/12
f[7] = -q*L/2
f[8] = q*L**2/12



print(K)
print(f"f={f}")


free_DOFs = [3,4,5,6,7,8,9,11]
constrained_DOFs = [0,1,2,10]

Kff = K[ix_(free_DOFs,free_DOFs)]
Kfc = K[ix_(free_DOFs,constrained_DOFs)]
Kcf = K[ix_(constrained_DOFs,free_DOFs)]
Kcc = K[ix_(constrained_DOFs,constrained_DOFs)]

ff= f[free_DOFs]
fc = f[constrained_DOFs]

u = zeros((NDOFs,1 ))

u[free_DOFs] = solve(Kff, ff)

R = Kcf@u[free_DOFs] + Kcc@u[constrained_DOFs] - fc

print(f"u={u}")
print(f"R={R}")
