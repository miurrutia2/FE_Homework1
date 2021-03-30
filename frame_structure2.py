from numpy import array, pi, zeros, ix_
from beam_element import beam_element
from scipy.linalg import solve
from math import sqrt


xy = array([
	[0,0],
	[0,3],
	[0,5],
	[0,6],
	[6,6.5],
	[6,5],
	[6,3],
	[6,0]
	])

conec = array([
	[0,1],
	[1,2],
	[2,3],
	[3,4],
	[4,5],
	[5,6],
	[6,7],
	[2,5],
	[1,6]
	], dtype = int)

peso_propio = -2400. #kg/m^3
properties_beam = {}
properties_roof = {}
properties_columns1 = {}
properties_columns2 = {}
properties_columns3 = {}
properties_columns5 = {}
properties_columns6 = {}
properties_columns7 = {}

properties_beam["E"] = 30e8
properties_beam["A"] = 0.2*0.4 #m^2
properties_beam["I"] = (0.2*(0.4**3))/12

properties_beam["qx"] = 0
properties_beam["qy"] = peso_propio*0.2*0.4*6

properties_columns1["E"] = 30e8
properties_columns1["A"] = 0.3*0.3 #m^2
properties_columns1["I"] = (0.3*(0.3**3))/12

properties_columns1["qx"] = peso_propio*0.3*0.3*3
properties_columns1["qy"] = 0

properties_columns2["E"] = 30e8
properties_columns2["A"] = 0.3*0.3 #m^2
properties_columns2["I"] = (0.3*(0.3**3))/12

properties_columns2["qx"] = peso_propio*0.3*0.3*2
properties_columns2["qy"] = 0

properties_columns3["E"] = 30e8
properties_columns3["A"] = 0.3*0.3 #m^2
properties_columns3["I"] = (0.3*(0.3**3))/12

properties_columns3["qx"] = peso_propio*0.3*0.3*1
properties_columns3["qy"] = 0

properties_columns5["E"] = 30e8
properties_columns5["A"] = 0.3*0.3 #m^2
properties_columns5["I"] = (0.3*(0.3**3))/12

properties_columns5["qx"] = peso_propio*0.3*0.3*1.5
properties_columns5["qy"] = 0

properties_columns7 = properties_columns1
properties_columns6 = properties_columns2

properties_roof["E"] = 30e8
properties_roof["A"] = 0.2*0.2 #m^2
properties_roof["I"] = (0.2*(0.2**3))/12

properties_roof["qx"] = 0
properties_roof["qy"] = peso_propio*0.2*0.2*6


properties = [properties_columns1,properties_columns2,properties_columns3,properties_roof,properties_columns5,properties_columns6,properties_columns7,properties_beam,properties_beam]

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

q = peso_propio
print(f"ke = {ke}")
L1 = 3.
L2 = 2.
L3 = 1.
L4 = sqrt(6**2 + 0.5**2)
L5 = 1.5
L6 = 2.
L7 = 3.
L8 = 6.
L9 = 6.

q1 = peso_propio * 0.3 * 0.3 * L1
q2 = peso_propio * 0.3 * 0.3 * L2
q3 = peso_propio * 0.3 * 0.3 * L3
q4 = peso_propio * 0.2 * 0.2 * L4
q5 = peso_propio * 0.3 * 0.3 * L5
q6 = peso_propio * 0.3 * 0.3 * L6
q7 = peso_propio * 0.3 * 0.3 * L7
q8 = peso_propio * 0.2 * 0.4 * L8
q9 = peso_propio * 0.2 * 0.4 * L9



f[0] = -q1*L1/2
f[1] = 0
f[2] = q1*L1**2/12
f[3] = -q1*L1/2 - q2*L2/2
f[4] = q9*L9/2
f[5] = -q1*L1**2/12 + q2*L2**2/12 + q9*L9**2/12
f[6] = -q2*L2/2 - q*L3/2
f[7] = q8*L8/2
f[8] = -q2*L2**2/12 + q3*L3**2/12 + q8*L8**2/12
f[9] = -q3*L3/2
f[10] = q4*L4/2
f[11] = -q3*L3**2/12 + q4*L4**2/12
f[12] = -q5*L5/2
f[13] = q4*L4/2
f[14] = -q5*L5**2/12 - q4*L4**2/12
f[15] = -q6*L6/2 - q5*L5/2
f[16] = q8*L8/2
f[17] = (-q6*L6**2 + q5*L5**2 - q8*L8**2)/12
f[18] = -q7*L7/2 - q6*L6/2
f[19] = q9*L9/2
f[20] = (-q7*L7**2 + q6*L6**2 - q9*L9**2)/12
f[21] = -q7*L7/2
f[22] = 0
f[23] = q7*L7**2/12


print(K)
print(f)


free_DOFs = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
constrained_DOFs = [0,1,2,21,22,23]

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