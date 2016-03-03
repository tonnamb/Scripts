import os
print os.getcwd()

import shutil, errno

import numpy as np
import math
import scipy as sp
from pylab import plot
from scipy.stats import norm

VaryTime = range(0,3001,1000)

RCBot = 5.0
RCTop = 32.0

runrange = range(1,51) + range(70,91)

def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise
		
if not os.path.exists("ui/data"):
	os.makedirs("ui/data")

for i in runrange:
    with open(str(i) + "_window/out.colvars.traj") as file:
        line = file.readline()
        print i
	line = file.readline()
        with open("ui/data/" + str(i),"w") as out:
            while line:
                rows = line.split()
		if rows[0] != '#':
			out.write(rows[1] + ' 299.79 ' + rows[2] + '\n')
		line = file.readline()

os.chdir("ui")		

# Check length
test = np.loadtxt("data/1")
datalength = len(test)

def factors(n):    
    return set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

# Initialization
Nwin = len(runrange)
Winkeep = runrange
data = np.empty([Nwin,datalength])
for i in range(Nwin):
	data[i] = np.genfromtxt("data/" + str(Winkeep[i]), usecols=(0))

# Definition of MK_Trend test
def mk_test(x, alpha = 0.05):
	n = len(x)
	# calculate S 
	s = 0
	for k in xrange(n-1):
		for j in xrange(k+1,n):
			s += np.sign(x[j] - x[k])
	
	# calculate the unique data
	unique_x = np.unique(x)
	g = len(unique_x)
	
	# calculate the var(s)
	if n == g: # there is no tie
		var_s = (n*(n-1)*(2*n+5))/18
	else: # there are some ties in data
		tp = np.zeros(unique_x.shape)
		for i in xrange(len(unique_x)):
			tp[i] = sum(unique_x[i] == x)
		var_s = (n*(n-1)*(2*n+5) + np.sum(tp*(tp-1)*(2*tp+5)))/18
	
	if s>0:
		z = (s - 1)/np.sqrt(var_s)
	elif s == 0:
			z = 0
	elif s<0:
		z = (s + 1)/np.sqrt(var_s)
	
	# calculate the p_value
	p = 2*(1-norm.cdf(abs(z))) # two tail test
	h = abs(z) > norm.ppf(1-alpha/2) 
	
	return h, p

# Definition of von Neumann test for Serial Correlation
def vn_test(x):
		sum = 0
		n = len(x)
		for i in range(n-1):
			sum = sum + (x[i+1]-x[i])**2
		qsquared = sum/(2*(n-1.0))
		r = qsquared/np.var(x)
		sigmar = ((1.0+1.0/(n-1.0))/(n+1.0))**0.5
		ur = (r-1.0)/sigmar
		p = 1.0-norm.cdf(abs(ur)) # 1-tail test
		return p

# Varying starting time of production

mkmean_t = np.empty(len(VaryTime))
mkstd_t = np.empty(len(VaryTime))
swtest_t = np.empty(len(VaryTime))
vntest_t = np.empty(len(VaryTime))
MinVal_t = np.empty(len(VaryTime))
MinNumSeg_t = np.empty(len(VaryTime))
MinWidth_t = np.empty(len(VaryTime))

for StartTimeIdx in range(len(VaryTime)):
	
	print VaryTime[StartTimeIdx]
	
	# Adjusting start of production
	DataStart = data[:,(VaryTime[StartTimeIdx]):]
	
	# MinSegWidth = int(math.ceil(len(DataStart[0])/100.0))
	# MaxSegWidth = int(math.floor(len(DataStart[0])/25.0))
	
	# Choosing # of segments
	lengthfactor = factors(len(DataStart[0]))
	NumSeg = sorted([y for y in lengthfactor if y>=30 and y<=100])
			
	# Test completion
			
	mkmean = np.empty([len(NumSeg),Nwin])
	mkstd = np.empty([len(NumSeg),Nwin])
	swtest = np.empty([len(NumSeg),Nwin])
	vntest = np.empty([len(NumSeg),Nwin])
	for seg in range(len(NumSeg)):
		testmean = np.mean(np.reshape(DataStart,(Nwin,NumSeg[seg],-1)),axis=-1)
		teststd = np.std(np.reshape(DataStart,(Nwin,NumSeg[seg],-1)),axis=-1)
		for i in range(len(testmean)):
			mkmean[seg,i] = mk_test(testmean[i],0.05)[1]
			mkstd[seg,i] = mk_test(teststd[i],0.05)[1]
			swtest[seg,i] = sp.stats.shapiro(testmean[i])[1]
			vntest[seg,i] = vn_test(testmean[i])
	
	mkmean_h = np.empty([len(NumSeg),Nwin])
	mkstd_h = np.empty([len(NumSeg),Nwin])
	swtest_h = np.empty([len(NumSeg),Nwin])
	vntest_h = np.empty([len(NumSeg),Nwin])
	
	for seg in range(len(NumSeg)):
		for i in range(len(testmean)):
			mkmean_h[seg,i] = int(mkmean[seg,i] < 0.05)
			mkstd_h[seg,i] = int(mkstd[seg,i] < 0.05)
			swtest_h[seg,i] = int(swtest[seg,i] < 0.05)
			vntest_h[seg,i] = int(vntest[seg,i] < 0.05)
			
	mkmean_s = np.sum(mkmean_h,axis=1)
	mkstd_s = np.sum(mkstd_h,axis=1)
	swtest_s = np.sum(swtest_h,axis=1)
	vntest_s = np.sum(vntest_h,axis=1)
	alltest = mkmean_s + mkstd_s + swtest_s + vntest_s
	
	val, idx = min((val, idx) for (idx, val) in enumerate(alltest))
	
	mkmean_t[StartTimeIdx] = mkmean_s[idx]
	mkstd_t[StartTimeIdx]  = mkstd_s[idx]
	swtest_t[StartTimeIdx] = swtest_s[idx]
	vntest_t[StartTimeIdx] = vntest_s[idx]
	MinVal_t[StartTimeIdx] = val
	MinNumSeg_t[StartTimeIdx] = NumSeg[idx]
	MinWidth_t[StartTimeIdx]  = len(DataStart[0])/NumSeg[idx]

val, idx = min((val, idx) for (idx, val) in enumerate(MinVal_t))
	
with open("stats.txt","w") as out:
	out.write("Use SegWidth -seg:\n")
	out.write(str(int(MinWidth_t[idx])) + "\n\n")
	
	out.write("Use SkipStart -ss:\n")
	out.write(str(VaryTime[idx]) + "\n\n")
	
	out.write("Best All Fail %:\n")
	out.write("{0:.1f}".format(100*val/(4*Nwin)) + "\n\n")	
	
	out.write("Best All Fail Count:\n")
	out.write(str(int(val)) + "\n\n")
	
	out.write("Best MK Mean Fail Count:\n")
	out.write(str(int(mkmean_t[idx])) + "\n\n")
	out.write("Best MK Std Fail Count:\n")
	out.write(str(int(mkstd_t[idx])) + "\n\n")
	out.write("Best SW Fail Count:\n")
	out.write(str(int(swtest_t[idx])) + "\n\n")
	out.write("Best VN Fail Count:\n")
	out.write(str(int(vntest_t[idx])) + "\n\n")
	
	out.write("Start Time Vary:\n")
	out.write(str(VaryTime) + "\n\n")
	
	out.write("Trend All Fail %:\n")
	out.write(str(np.round((100*MinVal_t/(4*Nwin)),1)) + "\n\n")	
	
os.system("ui.out -ui -T 300 -min 0 -max 50 -n 2000 -u kcal -seg " + str(int(MinWidth_t[idx])) + " -ss " + str(VaryTime[idx]) + " -r -1 -v 2 > log.txt")

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

found_min1 = find_nearest(rc_val,RCBot)
found_min2 = find_nearest(rc_val,RCTop)
id1 = rc_idx[found_min1[1]]
id2 = rc_idx[found_min2[1]]

filter1 = np.where(de_id2==id1)

id_dif = np.argwhere(de_id1[filter1[0]]==id2)[0,0]+filter1[0][0]


with open("flux_barrier.txt",'w') as out:
	out.write("dA = " + str(de_val[id_dif]) + " +- " + str(de_std[id_dif]) + " kcal/mol\n")
	out.write("From RC " + str(rc_val[found_min2[1]]) + " to " + str(rc_max_val[found_min1[1]]) + " (" + str(id2) + "-" + str(id1) + ")\n\n")