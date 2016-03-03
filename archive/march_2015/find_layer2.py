import os
print os.getcwd()

import numpy as np
from numpy import *
from numpy import linalg as LA

# Initialization

fname = "data.txt"
DataXYZ = genfromtxt(fname)

# Define new NAxis

b1 = -0.47
b2 = -0.55
b3 = -0.7

NAxis = np.array([b1,b2,b3])/LA.norm([b1,b2,b3])

# Project data to new NAxis

DataDot = dot(DataXYZ,NAxis)

# Plot histogram of layers
# import pylab
# pylab.hist(DataDot, bins = 200, normed =1)
# pylab.show()

# Locate atoms in layer

a1 = -45.0
a2 = -41.0
a3 = -38.5
a4 = -36.0
a5 = -32.5

layer1 = np.where(logical_and( DataDot < a2, DataDot > a1 ))
layer2 = np.where(logical_and( DataDot < a3, DataDot > a2 ))
layer3 = np.where(logical_and( DataDot < a4, DataDot > a3 ))
layer4 = np.where(logical_and( DataDot < a5, DataDot > a4 ))


# Adjust to LAMMPS numbering

layer1 = layer1[0] + 15811
layer2 = layer2[0] + 15811
layer3 = layer3[0] + 15811
layer4 = layer4[0] + 15811

# Create directory
if not os.path.exists("layer"):
	os.makedirs("layer")

# Write into list
	
f = file('layer/layer1.txt','w')
np.savetxt(f,layer1,fmt="%d",newline=' ')
f.close()

f = file('layer/layer2.txt','w')
np.savetxt(f,layer2,fmt="%d",newline=' ')
f.close()

f = file('layer/layer3.txt','w')
np.savetxt(f,layer3,fmt="%d",newline=' ')
f.close()

f = file('layer/layer4.txt','w')
np.savetxt(f,layer4,fmt="%d",newline=' ')
f.close()