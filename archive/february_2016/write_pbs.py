import os
import shutil

	
def writerun(run,start,stop,kforce):

	f = open("1ns/k_" + str(kforce) + "/" + str(run) + "pbs_run.txt","w")

	f.write("#PBS -l nodes=8\n")
	f.write("#PBS -l walltime=24:00:00\n")
	f.write("#PBS -j oe\n")
	f.write("#PBS -M tub179@psu.edu\n")
	f.write("#PBS -m bae\n")
	f.write("#PBS -N " + str(run) + "_synPVP10_100_1ns\n")
	f.write("#PBS -l pmem=2gb\n")
	f.write("#PBS -q lionxj-kaf2\n\n")

	f.write("cd $PBS_O_WORKDIR\n\n")

	f.write("echo \" \"\n")
	f.write("echo \"Starting job on `hostname` at `date`\"\n")
	f.write("echo \" \"\n\n")

	f.write("module load openmpi/intel/1.6.5\n")
	f.write("module load lammps/17Dec13\n")
	f.write("module load python\n\n")

	f.write("PATH=\"/gpfs/home/tub179/work/Research:$PATH\"\n\n")

	f.write("python run_window" + str(run) + ".py\n\n")

	f.write("echo \" \"\n")
	f.write("echo \"Completing job on `hostname` at `date`\"\n")
	f.write("echo \" \"\n")

	f.close()

for i in range(1,102):
	writerun(i,i,i,0.7)