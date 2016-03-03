import os
import shutil

def writesystem(window,windowstep,kforce):
    
	os.makedirs("2ns/k_" + str(kforce) + "/" + window + "_window")
    
	shutil.copy("Ag_O1X5.5_O2X0.55.eam.fs","2ns/k_" + str(kforce) + "/" + window + "_window" + "/Ag_O1X5.5_O2X0.55.eam.fs")
	
	f = open("2ns/k_" + str(kforce) + "/" + window + "_window/system.in","w")
	
	f.write("\n# ----------------- Settings Section -----------------\n\n")
	
	f.write("units           metal\n")
	f.write("boundary        p p p\n")
	f.write("neigh_modify    delay 0 every 1 check yes\n\n")
	
	f.write("atom_style      full\n")
	f.write("bond_style      harmonic\n")
	f.write("angle_style     charmm\n")
	f.write("dihedral_style  charmm\n")
	f.write("improper_style  harmonic\n")
	f.write("pair_style      hybrid/overlay eam/fs lj/charmm/coul/long 10.0 12.0 zhou 12.0 0.75 20.0 morse 5.5\n")
	f.write("kspace_style    pppm 1e-4\n\n")
	
	f.write("# ----------------- Atom Definition Section -----------------\n\n")
	
	f.write("read_restart US.window" + window + ".restart.1ns\n\n")
	
	f.write("pair_coeff   *   *   eam/fs Ag_O1X5.5_O2X0.55.eam.fs NULL NULL NULL NULL NULL O1 NULL NULL NULL O2 Ag\n")
	f.write("pair_coeff  1*2   11    morse   0.00480  1.30   1.94 # HGA1, HGA2, HGA3\n")
	f.write("pair_coeff    3   11    morse   0.00722  1.14   4.85 # CG2R53\n")
	f.write("pair_coeff    4   11    morse   0.00230  1.03   4.92 # CG3C52\n")
	f.write("pair_coeff    5   11    morse   0.00492  1.76   3.26 # NG2R53\n")
	f.write("pair_coeff    6   11    morse   0.00325  3.34   2.65 # OG2D1\n")
	f.write("pair_coeff    7   11    morse   0.00480  1.30   1.94 # HCA2\n")
	f.write("pair_coeff    8   11    morse   0.00448  1.06   2.13 # HCP1\n")
	f.write("pair_coeff    9   11    morse   0.00216  1.26   4.99 # CC32A\n")
	f.write("pair_coeff   10   11    morse   0.00217  2.09   3.43 # OC311\n")
	f.write("pair_coeff  1*2   11    zhou    0.0 1.0 1.0 10.2847 2.361 # H\n")
	f.write("pair_coeff  3*4   11    zhou    0.0 1.0 1.0 36.3619 2.812 # C\n")
	f.write("pair_coeff    5   11    zhou    0.0 1.0 1.0 30.4846 2.757 # N\n")
	f.write("pair_coeff    6   11    zhou    0.0 1.0 1.0 22.9973 2.702 # O\n")
	f.write("pair_coeff    7   11    zhou    0.0 1.0 1.0 10.2847 2.361 # HCA2\n")
	f.write("pair_coeff    8   11    zhou    0.0 1.0 1.0 10.2847 2.361 # HCP1\n")
	f.write("pair_coeff    9   11    zhou    0.0 1.0 1.0 36.3619 2.812 # CC32A\n")
	f.write("pair_coeff    10  11    zhou    0.0 1.0 1.0 22.9973 2.702 # OC311\n")
	f.write("pair_coeff   1   1  lj/charmm/coul/long    0.001518  2.387609  0.001518  2.387609 #HGA2   HGA2\n")
	f.write("pair_coeff   1   2  lj/charmm/coul/long    0.001740  1.393811  0.001740  1.393811 #HGA2   HGP1\n")
	f.write("pair_coeff   1   3  lj/charmm/coul/long    0.001147  3.153782  0.001147  3.153782 #HGA2   CG2R53\n")
	f.write("pair_coeff   1   4  lj/charmm/coul/long    0.001987  2.993420  0.000811  2.886512 #HGA2   CG3C52\n")
	f.write("pair_coeff   1   5  lj/charmm/coul/long    0.003628  2.841967  0.003628  2.841967 #HGA2   NG2R53\n")
	f.write("pair_coeff   1   6  lj/charmm/coul/long    0.002810  2.708333  0.002810  2.441063 #HGA2   OG2D1\n")
	f.write("pair_coeff   1   7  lj/charmm/coul/long    0.001518  2.387609  0.001518  2.387609 #HGA2   HCA2\n")
	f.write("pair_coeff   1   8  lj/charmm/coul/long    0.001740  1.393811  0.001740  1.393811 #HGA2   HCP1\n")
	f.write("pair_coeff   1   9  lj/charmm/coul/long    0.001920  2.984511  0.000811  2.886512 #HGA2   CC32A\n")
	f.write("pair_coeff   1  10  lj/charmm/coul/long    0.003556  2.766241  0.003556  2.766241 #HGA2   OC311\n")
	f.write("pair_coeff   2   2  lj/charmm/coul/long    0.001995  0.400014  0.001995  0.400014 #HGP1   HGP1\n")
	f.write("pair_coeff   2   3  lj/charmm/coul/long    0.001315  2.159984  0.001315  2.159984 #HGP1   CG2R53\n")
	f.write("pair_coeff   2   4  lj/charmm/coul/long    0.002278  1.999622  0.000930  1.892714 #HGP1   CG3C52\n")
	f.write("pair_coeff   2   5  lj/charmm/coul/long    0.004159  1.848169  0.004159  1.848169 #HGP1   NG2R53\n")
	f.write("pair_coeff   2   6  lj/charmm/coul/long    0.003222  1.714535  0.003222  1.447265 #HGP1   OG2D1\n")
	f.write("pair_coeff   2   7  lj/charmm/coul/long    0.001740  1.393811  0.001740  1.393811 #HGP1   HCA2\n")
	f.write("pair_coeff   2   8  lj/charmm/coul/long    0.001995  0.400014  0.001995  0.400014 #HGP1   HCP1\n")
	f.write("pair_coeff   2   9  lj/charmm/coul/long    0.002201  1.990713  0.000930  1.892714 #HGP1   CC32A\n")
	f.write("pair_coeff   2  10  lj/charmm/coul/long    0.004076  1.772443  0.004076  1.772443 #HGP1   OC311\n")
	f.write("pair_coeff   3   3  lj/charmm/coul/long    0.000867  3.919954  0.000867  3.919954 #CG2R53 CG2R53\n")
	f.write("pair_coeff   3   4  lj/charmm/coul/long    0.001502  3.759593  0.000613  3.652685 #CG2R53 CG3C52\n")
	f.write("pair_coeff   3   5  lj/charmm/coul/long    0.002743  3.608140  0.002743  3.608140 #CG2R53 NG2R53\n")
	f.write("pair_coeff   3   6  lj/charmm/coul/long    0.002124  3.474505  0.002124  3.207235 #CG2R53 OG2D1\n")
	f.write("pair_coeff   3   7  lj/charmm/coul/long    0.001147  3.153782  0.001147  3.153782 #CG2R53 HCA2\n")
	f.write("pair_coeff   3   8  lj/charmm/coul/long    0.001315  2.159984  0.001315  2.159984 #CG2R53 HCP1\n")
	f.write("pair_coeff   3   9  lj/charmm/coul/long    0.001451  3.750684  0.000613  3.652685 #CG2R53 CC32A\n")
	f.write("pair_coeff   3  10  lj/charmm/coul/long    0.002688  3.532413  0.002688  3.532413 #CG2R53 OC311\n")
	f.write("pair_coeff   4   4  lj/charmm/coul/long    0.002602  3.599231  0.000434  3.385415 #CG3C52 CG3C52\n")
	f.write("pair_coeff   4   5  lj/charmm/coul/long    0.004750  3.447778  0.001939  3.340870 #CG3C52 NG2R53\n")
	f.write("pair_coeff   4   6  lj/charmm/coul/long    0.003680  3.314144  0.001502  2.939966 #CG3C52 OG2D1\n")
	f.write("pair_coeff   4   7  lj/charmm/coul/long    0.001987  2.993420  0.000811  2.886512 #CG3C52 HCA2\n")
	f.write("pair_coeff   4   8  lj/charmm/coul/long    0.002278  1.999622  0.000930  1.892714 #CG3C52 HCP1\n")
	f.write("pair_coeff   4   9  lj/charmm/coul/long    0.002514  3.590322  0.000434  3.385415 #CG3C52 CC32A\n")
	f.write("pair_coeff   4  10  lj/charmm/coul/long    0.004655  3.372052  0.001901  3.265144 #CG3C52 OC311\n")
	f.write("pair_coeff   5   5  lj/charmm/coul/long    0.008673  3.296325  0.008673  3.296325 #NG2R53 NG2R53\n")
	f.write("pair_coeff   5   6  lj/charmm/coul/long    0.006718  3.162691  0.006718  2.895421 #NG2R53 OG2D1\n")
	f.write("pair_coeff   5   7  lj/charmm/coul/long    0.003628  2.841967  0.003628  2.841967 #NG2R53 HCA2\n")
	f.write("pair_coeff   5   8  lj/charmm/coul/long    0.004159  1.848169  0.004159  1.848169 #NG2R53 HCP1\n")
	f.write("pair_coeff   5   9  lj/charmm/coul/long    0.004589  3.438869  0.001939  3.340870 #NG2R53 CC32A\n")
	f.write("pair_coeff   5  10  lj/charmm/coul/long    0.008500  3.220599  0.008500  3.220599 #NG2R53 OC311\n")
	f.write("pair_coeff   6   6  lj/charmm/coul/long    0.005204  3.029056  0.005204  2.494516 #OG2D1  OG2D1\n")
	f.write("pair_coeff   6   7  lj/charmm/coul/long    0.002810  2.708333  0.002810  2.441063 #OG2D1  HCA2\n")
	f.write("pair_coeff   6   8  lj/charmm/coul/long    0.003222  1.714535  0.003222  1.447265 #OG2D1  HCP1\n")
	f.write("pair_coeff   6   9  lj/charmm/coul/long    0.003555  3.305235  0.001502  2.939966 #OG2D1  CC32A\n")
	f.write("pair_coeff   6  10  lj/charmm/coul/long    0.006584  3.086964  0.006584  2.819694 #OG2D1  OC311\n")
	f.write("pair_coeff   7   7  lj/charmm/coul/long    0.001518  2.387609  0.001518  2.387609 #HCA2   HCA2\n")
	f.write("pair_coeff   7   8  lj/charmm/coul/long    0.001740  1.393811  0.001740  1.393811 #HCA2   HCP1\n")
	f.write("pair_coeff   7   9  lj/charmm/coul/long    0.001920  2.984511  0.000811  2.886512 #HCA2   CC32A\n")
	f.write("pair_coeff   7  10  lj/charmm/coul/long    0.003556  2.766241  0.003556  2.766241 #HCA2   OC311\n")
	f.write("pair_coeff   8   8  lj/charmm/coul/long    0.001995  0.400014  0.001995  0.400014 #HCP1   HCP1\n")
	f.write("pair_coeff   8   9  lj/charmm/coul/long    0.002201  1.990713  0.000930  1.892714 #HCP1   CC32A\n")
	f.write("pair_coeff   8  10  lj/charmm/coul/long    0.004076  1.772443  0.004076  1.772443 #HCP1   OC311\n")
	f.write("pair_coeff   9   9  lj/charmm/coul/long    0.002428  3.581413  0.000434  3.385415 #CC32A  CC32A\n")
	f.write("pair_coeff   9  10  lj/charmm/coul/long    0.004498  3.363143  0.001901  3.265144 #CC32A  OC311\n")
	f.write("pair_coeff  10  10  lj/charmm/coul/long    0.008330  3.144872  0.008330  3.144872 #OC311  OC311 \n\n")
	
	f.write("special_bonds   charmm\n\n")
	
	f.write("# ----------------- Run Section -----------------\n\n")
	
	f.write("timestep        0.001\n\n")
	
	f.write("restart         2000 US.window" + window + ".restart.2ns_dup US.window" + window + ".restart.2ns\n\n")
	
	f.write("thermo          0\n")
	f.write("thermo_style    custom step pe ke etotal temp vol press lx lz pxx pyy pzz\n")
	f.write("thermo_modify   flush yes\n\n")
	
	f.write("variable STEP equal step\n")
	f.write("variable PE equal pe\n")
	f.write("variable KE equal ke\n")
	f.write("variable ETOTAL equal etotal\n")
	f.write("variable TEMP equal temp\n")
	f.write("variable VOL equal vol\n")
	f.write("variable PRESS equal press\n")
	f.write("variable LX equal lx\n")
	f.write("variable LZ equal lz\n")
	f.write("variable PXX equal pxx\n")
	f.write("variable PYY equal pyy\n")
	f.write("variable PZZ equal pzz\n\n")
	
	f.write("fix            1 all nvt temp 300.0 300.0 0.1 fixedpoint 0.0 0.0 0.0\n")
	f.write("fix			2 all colvars pull_unix.in\n")
	f.write("fix 			thermo_output all print 100 \"${STEP} ${PE} ${KE} ${ETOTAL} ${TEMP} ${VOL} ${PRESS} ${LX} ${LZ} ${PXX} ${PYY} ${PZZ}\" file thermo.lammps title \"# Step PotEng KinEng TotEng Temp Volume Press Lx Lz Pxx Pyy Pzz\" screen no\n\n")
	
	f.write("dump		    1 all custom 20000 traj_npt.lammpstrj id mol type x y z ix iy iz\n")
	
	f.write("run             1000000\n")
	f.close()
	
	f = open("2ns/k_" + str(kforce) + "/" + window + "_window/pull_unix.in","w")
	
	f.write("colvarsTrajFrequency 100")
	f.write("\n")
	f.write("colvar {\n")
	f.write("	name ZPULL\n")
	f.write("	width 1.0\n")
	f.write("	UpperBoundary 46.0\n")
	f.write("	LowerBoundary 0.0\n")
	f.write("	distanceZ {\n")
	f.write("		main { atomNumbersRange 3433-3445 }\n")
	f.write("		ref { \n")
	f.write("			atomNumbersRange 3001-3036\n")
	f.write("			atomNumbersRange 3109-3144\n")
	f.write("		}\n")
	f.write("		axis (0.0,0.0,-1.0)\n")
	f.write("	}\n")
	f.write("}\n")
	f.write("\n")
	f.write("harmonic {\n")
	f.write("	colvars ZPULL\n")
	f.write("	forceConstant " + str(kforce) + "\n")    
	f.write("	centers " + str((float(window) - 1.0)*windowstep+2) + "\n")
	f.write("	targetCenters " + str((float(window) - 1.0)*windowstep+2) + "\n")
	f.write("	targetNumSteps 1000000\n")
	f.write("	outputCenters yes\n")
	f.write("	outputAccumulatedWork yes\n")
	f.write("}\n")
	
	f.close()
	
	f = open("2ns/k_" + str(kforce) + "/" + window + "_window/run_npt.py","w")
	
	f.write("import os\n")
	f.write("os.system(\"mpirun lmp_openmpi_17Dec13_Grimme_s_d < system.in\")\n")
	
	f.close()	
	
	shutil.copy("1ns/k_" + str(kforce) + "/" + window + "_window/US.window" + window + ".restart.1ns","2ns/k_" + str(kforce) + "/" + window + "_window/US.window" + window + ".restart.1ns")


	
def writerun(run,start,stop,kforce):

	f = open("2ns/k_" + str(kforce) + "/" + str(run) + "pbs_run.txt","w")

	f.write("#PBS -l nodes=8\n")
	f.write("#PBS -l walltime=24:00:00\n")
	f.write("#PBS -j oe\n")
	f.write("#PBS -M tub179@psu.edu\n")
	f.write("#PBS -m bae\n")
	f.write("#PBS -N " + str(run) + "_k" + str(kforce) + "_2ns_2P100_bs_US\n")
	f.write("#PBS -l pmem=2gb\n\n")

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

	f = open("2ns/k_" + str(kforce) + "/run_window" + str(run) + ".py","w")

	f.write("import os\n")
	f.write("import shutil\n")

	f.write("def runsystem(window):\n")
	f.write("	os.chdir(window + \"_window\")\n")
	f.write("	os.system(\"pwd\")\n")
	f.write("	os.system(\"mpirun lmp_openmpi_17Dec13_Grimme_s_d < system.in\")\n")
	f.write("	os.chdir(os.pardir)\n")

	f.write("for i in range(" + str(start) + "," + str(stop+1) + "):\n")
	f.write("	runsystem(str(i))\n")

	f.close()
	
for i in range(1,102):
	writesystem(str(i),0.18,0.7)

for i in range(1,34):
	writerun(i,i*3-2,i*3,0.7)

writerun(34,100,101,0.7)
