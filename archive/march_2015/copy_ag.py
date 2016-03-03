import os
import shutil

def writesystem(window):
    
	shutil.copy("Ag_O1X5.5_O2X0.55.eam.fs","1ns/" + window + "_window/Ag_O1X5.5_O2X0.55.eam.fs")
	
for i in range(1,102):
	writesystem(str(i))
