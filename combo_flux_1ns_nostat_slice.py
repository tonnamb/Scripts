import os
print os.getcwd()

import shutil, errno

import numpy as np
import math
import scipy as sp
from pylab import plot
from scipy.stats import norm
from itertools import islice

RCMax = 5.0
RCMin = 8.0

def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise
		
if not os.path.exists("ui_nostat/data_unslice"):
	os.makedirs("ui_nostat/data_unslice")

runrange = range(1,61)
    
for i in runrange:
    with open(str(i) + "_window/out.colvars.traj") as file:
        line = file.readline()
        print i
	line = file.readline()
        with open("ui_nostat/data_unslice/" + str(i),"w") as out:
            while line:
                rows = line.split()
		if rows[0] != '#':
			out.write(rows[1] + ' 39.20 ' + rows[2] + '\n')
		line = file.readline()

os.chdir("ui_nostat")

if not os.path.exists("data"):
	os.makedirs("data")

for i in runrange:
	with open('data_unslice/' + str(i)) as fin, open('data/' + str(i), 'w') as fout:
		fout.writelines(islice(fin, None, None, 10))


os.system("ui.out -ui -T 433 -min 0 -max 50 -n 2000 -u kcal -ss 1000 -r -1 -v 2 > log.txt")

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return np.array([array[idx],idx])

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

os.system("wham.out -wham -T 433 -min 0 -max 50 -n 400 -u kcal -ss 1000 -r -1 -v 2 > wham.txt")