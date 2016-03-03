import os
print os.getcwd()

import numpy as np
from numpy import *

# Initialization
fname = "comz.lammps"
DataComz = genfromtxt(fname, skip_header=1)
lz = 62.029937026475

# Slice DataComz into Size of Trajectory
DataComz = DataComz[20:-1:50]

# Wrap DataComz
DataComz[DataComz < 0] += lz
DataComz[DataComz > lz] -= lz

# Histogram to locate layers
# import pylab
# pylab.hist(DataComz[0], bins = 200, normed =1)
# pylab.show()
# Bottom layers:
# 1st layer = 21 - 25
# 2nd layer = 17 - 21
# 3rd layer = 13 - 17

# Locate firstlayer
firstlayer = []
for i in range(0,len(DataComz)):
	firstlayer.append(np.where(logical_and( DataComz[i] < 25, DataComz[i] > 21 )))

# Locate secondlayer
secondlayer = []
for i in range(0,len(DataComz)):
	secondlayer.append(np.where(logical_and( DataComz[i] < 21, DataComz[i] > 17 )))
	
# Locate thirdlayer
thirdlayer = []
for i in range(0,len(DataComz)):
	thirdlayer.append(np.where(logical_and( DataComz[i] < 17, DataComz[i] > 13 )))

# Create directory
if not os.path.exists("layer"):
	os.makedirs("layer")
	
# Save to text file
f = file('layer/firstlayer.txt','a')
for item in firstlayer:
	np.savetxt(f,item,fmt="%d")
f.close()

f = file('layer/secondlayer.txt','a')
for item in secondlayer:
	np.savetxt(f,item,fmt="%d")
f.close()

f = file('layer/thirdlayer.txt','a')
for item in thirdlayer:
	np.savetxt(f,item,fmt="%d")
f.close()

# Read text file to list
with open("layer/firstlayer.txt") as file:
	line = file.readline()
	with open("layer/serialfirstlayer.txt","w") as out:
		while line:
			rows = line.split()
			for i in rows:
				out.write(' '.join(str(e) for e in range(int(i)*10+1,int(i)*10+11)) + ' ')
			out.write('\n')
			line = file.readline()
			
with open("layer/secondlayer.txt") as file:
	line = file.readline()
	with open("layer/serialsecondlayer.txt","w") as out:
		while line:
			rows = line.split()
			for i in rows:
				out.write(' '.join(str(e) for e in range(int(i)*10+1,int(i)*10+11)) + ' ')
			out.write('\n')
			line = file.readline()
			
with open("layer/thirdlayer.txt") as file:
	line = file.readline()
	with open("layer/serialthirdlayer.txt","w") as out:
		while line:
			rows = line.split()
			for i in rows:
				out.write(' '.join(str(e) for e in range(int(i)*10+1,int(i)*10+11)) + ' ')
			out.write('\n')
			line = file.readline()