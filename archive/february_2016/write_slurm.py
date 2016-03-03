import os
import shutil
	
def writerun(run,start,stop):

	f = open(str(run) + "slurm.txt","w")

	f.write("#!/bin/bash\n")
	f.write("#SBATCH -J " + str(run) + "_SDA_PVP10_100_1ns                       \n")
	f.write("#SBATCH -o " + str(run) + "_SDA_PVP10_100_1ns.o%j                   \n")
	f.write("#SBATCH -n 32                                      \n")
	f.write("#SBATCH -p normal                                 # queue (partition) -- normal, development, etc.\n")
	f.write("#SBATCH -t 24:00:00                                # run time (hh:mm:ss)\n")
	f.write("#SBATCH --mail-user=tub179@psu.edu\n")
	f.write("#SBATCH --mail-type=begin                         \n")
	f.write("#SBATCH --mail-type=end                           \n")
	f.write("#SBATCH -A TG-DMR110061\n")
	f.write("PATH=\"/home1/03207/tub179/exe:$PATH\"\n")
	f.write("python run_window" + str(run) + ".py\n\n")

	f.close()
	
	f = open("run_window" + str(run) + ".py","w")

	f.write("import os\n")
	f.write("import shutil\n")

	#f.write("def runsystem(window):\n")
	#f.write("	os.chdir(window + \"_window\")\n")
	#f.write("	os.system(\"pwd\")\n")
	#f.write("	os.system(\"ibrun lmp_stampede_omp < system.in\")\n")
	#f.write("	os.chdir(os.pardir)\n")
	
	f.write("def runsystem(window):\n")
	f.write("	os.chdir(window + \"_window/no_colvars\")\n")
	f.write("	os.system(\"ibrun lmp_stampede_omp < system.in\")\n")
	f.write("	os.chdir(os.pardir)\n")
	f.write("	shutil.copy(\"no_colvars/US.window\" + window + \".1ns_nc\",\"US.window\" + window + \".1ns_nc\")\n")
	f.write("	os.system(\"pwd\")\n")
	f.write("	os.system(\"ibrun lmp_stampede_omp < system.in\")\n")
	f.write("	os.chdir(os.pardir)\n")

	f.write("for i in range(" + str(start) + "," + str(stop+1) + "):\n")
	f.write("	runsystem(str(i))\n")

	f.close()

for i in range(1,102):
	writerun(i,i,i)
