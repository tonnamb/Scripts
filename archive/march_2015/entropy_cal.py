import os
print os.getcwd()
import numpy as np
from scipy.integrate import simps

rc1 = 3.11369
loc1 = 31
rc2 = 28.2170150661
loc2 = 127

pe = []
for window in range(1,102):
	pe.append([])
	
print "Reading thermo"

nslist = ["1ns","2ns", "4ns", "6ns"]
for ns in nslist:
	print ns
	for i in range(1,102):
		with open(ns + "/k_0.7/" + str(i) + "_window/thermo.lammps") as file:
			print i
			line = file.readline()
			line = file.readline()
			while line:
				rows = line.split()
				pe[i-1].append(float(rows[1]))
				line = file.readline()
			print "window: " + str(i) + ", length: " + str(len(pe[i-1]))

print "Reading data for xval"	
				
xval = []
for window in range(1,102):
	xval.append([])
				
for i in range(1,102):
    with open(nslist[-1] + "/k_0.7/ui/data/" + str(i)) as file:
		print i
		line = file.readline()
		while line:
			rows = line.split()
			xval[i-1].append(float(rows[0]))
			line = file.readline()
			
rcval = []
dadxi = []
with open(nslist[-1] + "/k_0.7/ui/dadxi_norm.xy") as file:
	line = file.readline()
	line = file.readline()
	line = file.readline()
	while line:
		rows = line.split()
		rcval.append(float(rows[0]))
		dadxi.append(float(rows[1]))
		line = file.readline()
			

perc1 = []
perc2 = []	

#halfbinwidth = 0.0125 #0.5*50.0/2000.0
halfbinwidth = 0.13125 #0.5*105.0/400.0

rc1top = rc1 + halfbinwidth
rc1bot = rc1 - halfbinwidth
rc2top = rc2 + halfbinwidth
rc2bot = rc2 - halfbinwidth

print "Finding PotEng that is in the bins"

for i in range(1,102):
	print "window: " + str(i) + ", xval: " + str(len(xval[i-1])) + ", pe: " + str(len(pe[i-1]))
	for j in range(0,len(xval[i-1])):
		if rc1bot < xval[i-1][j] < rc1top:
			perc1.append(pe[i-1][j])
		if rc2bot < xval[i-1][j] < rc2top:
			perc2.append(pe[i-1][j])

print "For bound state: RC1 = " + str(rc1)
print "Number of observations: " + str(len(perc1))
print "Average PotE: " + str(sum(perc1)/len(perc1))
print "Max PotE: " + str(max(perc1))
print "Min PotE: " + str(min(perc1))
print "StDev PotE: " + str(np.std(perc1))

print "For free state: RC2 = " + str(rc2)
print "Number of observations: " + str(len(perc2))
print "Average PotE: " + str(sum(perc2)/len(perc2))
print "Max PotE: " + str(max(perc2))
print "Min PotE: " + str(min(perc2))
print "StDev PotE: " + str(np.std(perc2))	

with open("PotE_RC1.txt","w") as out:
	for z in perc1:
		out.write(str(z) + '\n')
		
with open("PotE_RC2.txt","w") as out:
	for z in perc2:
		out.write(str(z) + '\n')
		
print "Unified approach"

perc = []
for k in range(0,loc2-loc1):
	perc.append([])

for i in range(1,102):
	print "Window " + str(i)
	for j in range(0,len(xval[i-1])):
		for k in range(0,loc2-loc1):
			if (rcval[loc1+k-1]-halfbinwidth) < xval[i-1][j] < (rcval[loc1+k-1]+halfbinwidth):
				perc[k].append(pe[i-1][j])
				break

temp = 300				
boltz = 0.0019872041
eVtokcal = 23.0605
beta = 1/(temp*boltz)

#for k in range(0,len(perc)):
#	with open(str(k) + "_PotE.txt","w") as out:
#		for z in perc[k]:
#			out.write(str(z) + '\n')

ui = []
for i in range(0,len(perc)):
	ui.append(eVtokcal*sum(perc[i])/len(perc[i]))
	
dudxi = np.gradient(ui,2*halfbinwidth)

integrand = []
for i in range(0,len(perc)):
	integrand.append(dudxi[i]-ui[i]*dudxi[i])

#integrand1 = []
#integrand2 = []
#for i in range(0,loc2-loc1):
#	integrand1.append(dadxi[loc1+i-1]-(23.0605*sum(perc[i])/len(perc[i]))*dadxi[loc1+i-1])
#	integrand2.append(dadxi[loc1+i-1]+(23.0605*sum(perc[i])/len(perc[i]))*dadxi[loc1+i-1])

#dentropy1 = -simps(integrand1,even='first')/(boltz*temp)
#dentropy2 = -simps(integrand2,even='first')/(boltz*temp)

uexpsig = []
expsig = []
for z in perc[0]:
	uexpsig.append(eVtokcal*z*exp(-beta*(eVtokcal*z-ui[-1])))
	expsig.append(exp(-beta*(eVtokcal*z-ui[-1])))

du1 = np.mean(uexpsig)/np.mean(expsig)+ui[-1]
du2 = np.mean(uexpsig)/np.mean(expsig)-ui[-1]
	
dentropy = -simps(integrand,dx=halfbinwidth*2,even='first')/(boltz*temp)

with open("dS.txt","w") as out:
	out.write(str(dentropy) + '\n')