import os
print os.getcwd()

import numpy as np
import time

def factors(n):    
    return set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

SSVary = range(20000,40001,500)
RCNear = 32.0

for j in SSVary:

	os.makedirs("data_" + str(j))
	print j
	for i in range(1,102):
		with open("data/" + str(i),'r') as f:
			with open("data_" + str(j) + "/" + str(i),'w') as out:
				count = 0
				for line in f:
					if count > j:
						out.write(line)
					count = count + 1
					if count > j+20000-1:
						break

# Check length
# test = np.loadtxt("data/1")
# datalength = len(test)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return np.array([array[idx],idx])

for i in SSVary:

	print i
	
	# lengthfactor = factors(datalength-i)
	# NumSeg = sorted([y for y in lengthfactor if y>=40 and y<=100])
	# seg = (datalength-i)/NumSeg[0]
	# 
	# print seg

	start = time.time()
	os.system("ui.out -ui -d data_" + str(i) + " -T 300 -min -5 -max 100 -n 400 -u kcal -seg 625 -r -1 -v 0 > log_" + str(i) + ".txt")
	end = time.time()
	print end - start
	
	rc_idx = []
	rc_val = []
	de_id1 = []
	de_id2 = []
	de_val = []
	de_std = []

	with open('log_' + str(i) + '.txt') as f:
		for line in f:
			if len(line) == 26 and line[0:3] == 'Min':
				rc_idx.append(int(line[9:12]))
				rc_val.append(float(line[17:26]))
			if len(line) == 82:
				de_id1.append(int(line[27:30]))
				de_id2.append(int(line[32:35]))
				de_val.append(float(line[39:53]))
				de_std.append(float(line[60:73]))
				
	rc_idx = np.array(rc_idx)
	rc_val = np.array(rc_val)
	de_id1 = np.array(de_id1)
	de_id2 = np.array(de_id2)
	de_val = np.array(de_val)
	de_std = np.array(de_std)

	near = find_nearest(rc_val,RCNear)
	id1 = rc_idx[0]
	id2 = rc_idx[near[1]]

	filter1 = np.argwhere(de_id2==id1)
	id_dif = np.argwhere(de_id1[filter1]==id2)[0,0]
	
	with open("dA.txt",'a') as out:
		out.write("SS: " + str(i) + "\n")
		out.write("dA = " + str(de_val[id_dif]) + " +- " + str(de_std[id_dif]) + " kcal/mol\n")
		out.write("From RC " + str(rc_val[near[1]]) + " to " + str(rc_val[0]) + " (" + str(id2) + "-" + str(id1) + ")\n\n")
		
	with open("plotdA.txt",'a') as out:
		out.write(str(de_val[id_dif]) + "\n")
		
	with open("plotStd.txt",'a') as out:
		out.write(str(de_std[id_dif]) + "\n")
		
	with open("plotSS.txt",'a') as out:
		out.write(str(i) + "\n")
