import os
print os.getcwd()

import numpy as np

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return np.array([array[idx],idx])

RCMax = 6.0
RCMin = 8.0

rc_idx = []
rc_val = []
rc_max_idx = []
rc_max_val = []
de_id1 = []
de_id2 = []
de_val = []
de_std = []

with open('log.txt') as f:
	for line in f:
		if len(line) == 26 and line[0:3] == 'Min':
			rc_idx.append(int(line[9:12]))
			rc_val.append(float(line[17:26]))
		if len(line) == 26 and line[0:3] == 'Max':
			rc_max_idx.append(int(line[9:12]))
			rc_max_val.append(float(line[17:26]))
		if len(line) == 82:
			de_id1.append(int(line[27:30]))
			de_id2.append(int(line[32:35]))
			de_val.append(float(line[39:53]))
			de_std.append(float(line[60:73]))
			
rc_idx = np.array(rc_idx)
rc_val = np.array(rc_val)
rc_max_idx = np.array(rc_max_idx)
rc_max_val = np.array(rc_max_val)
de_id1 = np.array(de_id1)
de_id2 = np.array(de_id2)
de_val = np.array(de_val)
de_std = np.array(de_std)

found_min = find_nearest(rc_val,RCMin)
found_max = find_nearest(rc_max_val,RCMax)
id1 = rc_max_idx[found_max[1]]
id2 = rc_idx[found_min[1]]

filter1 = np.where(de_id2==id1)

id_dif = np.argwhere(de_id1[filter1[0]]==id2)[0,0]+filter1[0][0]


with open("flux_barrier.txt",'w') as out:
	out.write("dA = " + str(de_val[id_dif]) + " +- " + str(de_std[id_dif]) + " kcal/mol\n")
	out.write("From RC " + str(rc_val[found_min[1]]) + " to " + str(rc_max_val[found_max[1]]) + " (" + str(id2) + "-" + str(id1) + ")\n\n")