import os
print os.getcwd()

import numpy as np
import math
import scipy as sp
from pylab import plot
from scipy.stats import norm

# Check length
test = np.loadtxt("data/1")
datalength = len(test)

VaryTime = range(30000,50001,1000)

def factors(n):    
    return set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

# Initialization
Nwin = 61
Winkeep = range(1,21)+range(61,102)
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