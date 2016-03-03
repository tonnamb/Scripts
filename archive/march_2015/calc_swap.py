import os
print os.getcwd()

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("log.lammps",skip_header=3, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32))

swap = np.empty(len(data)-1)

for i in range(1,len(data)):
	count = 0
	for j in range(8):
		if data[i,j] != data[i-1,j]:
			count = count + 1
	swap[i-1] = count/2
	
maxswap = 4*(len(data)-1)
swaprate = np.sum(swap)/maxswap

plt.plot(range(1,len(swap)+1),swap)
plt.savefig('profile_swap.png')

with open("system.in","r") as fp:
	for i, line in enumerate(fp):
		if i == 166:
			temperature = line

with open("rate_swap.txt","w") as out:
	out.write("Temperature set:\n")
	out.write(temperature + "\n\n")
	
	out.write("Swap rate:\n")
	out.write(str(swaprate) + "\n\n")