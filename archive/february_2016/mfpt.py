import os
print os.getcwd()

import numpy as np
from scipy.integrate import simps

import matplotlib.pyplot as plt

pmf = []

with open('fe_ui.xy') as f:
	line = f.readline()
	line = f.readline()
	line = f.readline()
	while line:
		rows = line.split()
		rows = np.array(rows)		
		rows = rows.astype(np.float)
		pmf.append(rows)
		line = f.readline()

pmf = np.array(pmf)

d = 2.35E11 # Ang^2/s, Diffusion coeff of Ag in EG
x0 = pmf[0,0]
xi = pmf[np.argmin(pmf, axis=0)[1],0]
xipos = np.argmin(pmf, axis=0)[1]

xf = pmf[np.abs(pmf[:,0]-32).argmin(),0]
xfpos = np.abs(pmf[:,0]-32).argmin()

#xf = pmf[np.argmax(pmf, axis=0)[1],0]
#xfpos = np.argmax(pmf, axis=0)[1]

print "d = {0} Ang^2/s".format(d)
print "x0 = {0} Ang".format(x0)
print "xi = {0} Ang".format(xi)
print "xf = {0} Ang".format(xf)

# MFPT at x = xf
def f1(x):
	return np.exp(-x)

int1 = f1(pmf[0:xfpos,1])
int1num = simps(int1,pmf[0:xfpos,0])

def f2(x):
	return np.exp(x)/d

int2 = f2(pmf[xipos:xfpos,1])*int1num
mfpt_xf = simps(int2,pmf[xipos:xfpos,0])

print "MFPT at xf = {0}".format(mfpt_xf)

# MFPT at all x

mfpt = []

for i in range(xipos+1,xfpos+1):
	int1 = f1(pmf[0:i,1])
	int1num = simps(int1,pmf[0:i,0])
	int2 = f2(pmf[xipos:i,1])*int1num
	mfpt_i = simps(int2,pmf[xipos:i,0])
	mfpt.append(mfpt_i)

plt.plot(pmf[xipos:xfpos,0],mfpt)
plt.xlabel('z from surface (Ang)')
plt.ylabel('MFPT (second)')
plt.savefig('profile_mfpt.png')
