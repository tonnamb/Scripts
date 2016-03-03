import os
print os.getcwd()

import numpy as np
from numpy import *

# Initialization of data

atoms = []
bonds = []
angles = []
dihedrals = []
impropers = []
with open('anti_10mer_HRA2.data','r') as file:
	line = file.readline()
	while line:
		if 'Atoms' in line:
			break
		line = file.readline()
	line = file.readline()
	line = file.readline()
	while line:
		if line == '\n':
			break
		atoms.append(line.split())
		line = file.readline()
	line = file.readline()
	while line:
		if 'Bonds' in line:
			break
		line = file.readline()
	line = file.readline()
	line = file.readline()
	while line:
		if line == '\n':
			break
		bonds.append(line.split())
		line = file.readline()
	line = file.readline()
	while line:
		if 'Angles' in line:
			break
		line = file.readline()
	line = file.readline()
	line = file.readline()
	while line:
		if line == '\n':
			break
		angles.append(line.split())
		line = file.readline()
	line = file.readline()
	while line:
		if 'Dihedrals' in line:
			break
		line = file.readline()
	line = file.readline()
	line = file.readline()
	while line:
		if line == '\n':
			break
		dihedrals.append(line.split())
		line = file.readline()
	line = file.readline()
	while line:
		if 'Impropers' in line:
			break
		line = file.readline()
	line = file.readline()
	line = file.readline()
	while line:
		if line == '\n':
			break
		impropers.append(line.split())
		line = file.readline()
		
# Write .lt file
		
with open('lt_composed.txt','w') as out:
	out.write('PVP10 {\n\n')
	out.write('   write_once("Data Masses") {\n')
	out.write('\t@atom:1      1.008  # HRA2\n')
	out.write('\t@atom:2      1.008  # HGA1\n')
	out.write('\t@atom:3      1.008  # HGA2\n')
	out.write('\t@atom:4      1.008  # HGA3\n')
	out.write('\t@atom:5     12.011  # CG2R53\n')
	out.write('\t@atom:6     12.011  # CG311\n')
	out.write('\t@atom:7     12.011  # CG321\n')
	out.write('\t@atom:8     12.011  # CG331\n')
	out.write('\t@atom:9     12.011  # CG3C52\n')
	out.write('\t@atom:10    14.007  # NG2R53\n')
	out.write('\t@atom:11   15.9994  # OG2D1\n')
	out.write('  } \n\n')
	out.write('   write_once("Data Bond Coeffs") {\n')
	out.write('\t@bond:1   13.00902   1.530   # CG2R53   CG3C52 \n')
	out.write('\t@bond:2   19.94716   1.380   # CG2R53   NG2R53 \n')
	out.write('\t@bond:3   24.71714   1.235   # CG2R53   OG2D1  \n')
	out.write('\t@bond:4   9.64836    1.538   # CG311    CG321  \n')
	out.write('\t@bond:5   9.64836    1.538   # CG311    CG331  \n')
	out.write('\t@bond:6   13.39929   1.111   # CG311    HGA1 \n')
	out.write('\t@bond:7   13.87629   1.430   # CG311    NG2R53 \n')
	out.write('\t@bond:8   13.39929   1.111   # CG321    HGA2 \n')
	out.write('\t@bond:9   13.96301   1.111   # CG331    HGA3 \n')
	out.write('\t@bond:10  8.45586    1.530   # CG3C52   CG3C52 \n')
	out.write('\t@bond:11  13.31256   1.100   # CG3C52   HRA2 \n')
	out.write('\t@bond:12  16.04446   1.450   # CG3C52   NG2R53 \n')
	out.write('  }\n\n')
	out.write('   write_once("Data Angle Coeffs") {\n')
	out.write('\t@angle:1   3.03544	106.5  0.00000	0.000    # CG2R53   CG3C52   CG3C52\n')
	out.write('\t@angle:2   2.51508	111.0  0.00000	0.000    # CG2R53   CG3C52   HRA2\n')
	out.write('\t@angle:3   2.16817	120.0  0.00000	0.000    # CG2R53   NG2R53   CG311 \n')
	out.write('\t@angle:4   3.25226	111.0  0.00000	0.000    # CG2R53   NG2R53   CG3C52\n')
	out.write('\t@angle:5   2.53025	113.5  0.48394	2.561    # CG311    CG321    CG311 \n')
	out.write('\t@angle:6   1.44964	110.1  0.97698	2.179    # CG311    CG321    HGA2\n')
	out.write('\t@angle:7   1.44964	110.1  0.97698	2.179    # CG311    CG331    HGA3\n')
	out.write('\t@angle:8   2.16817	116.0  0.00000	0.000    # CG311    NG2R53   CG3C52\n')
	out.write('\t@angle:9   2.53025	113.5  0.48394	2.561    # CG321    CG311    CG321 \n')
	out.write('\t@angle:10  2.31344	114.0  0.34691	2.561    # CG321    CG311    CG331 \n')
	out.write('\t@angle:11  1.49604	110.1  0.97698	2.179    # CG321    CG311    HGA1\n')
	out.write('\t@angle:12  3.03544	113.5  0.00000	0.000    # CG321    CG311    NG2R53\n')
	out.write('\t@angle:13  1.49604	110.1  0.97698	2.179    # CG331    CG311    HGA1\n')
	out.write('\t@angle:14  3.03544	113.5  0.00000	0.000    # CG331    CG311    NG2R53\n')
	out.write('\t@angle:15  5.20361	105.5  0.00000	0.000    # CG3C52   CG2R53   NG2R53\n')
	out.write('\t@angle:16  2.81862	126.7  0.00000	0.000    # CG3C52   CG2R53   OG2D1 \n')
	out.write('\t@angle:17  2.51508	109.5  0.48394	2.561    # CG3C52   CG3C52   CG3C52\n')
	out.write('\t@angle:18  1.51772	111.4  0.97698	2.179    # CG3C52   CG3C52   HRA2\n')
	out.write('\t@angle:19  3.90271	104.5  0.00000	0.000    # CG3C52   CG3C52   NG2R53\n')
	out.write('\t@angle:20  2.08144	108.0  0.00000	0.000    # HGA1 CG311    NG2R53  \n')
	out.write('\t@angle:21  1.53940	109.0  0.23416	1.802    # HGA2 CG321    HGA2\n')
	out.write('\t@angle:22  1.53940	108.4  0.23416	1.802    # HGA3 CG331    HGA3\n')
	out.write('\t@angle:23  1.66949	106.8  0.23416	1.802    # HRA2 CG3C52   HRA2\n')
	out.write('\t@angle:24  2.55844	111.0  0.00000	0.000    # HRA2 CG3C52   NG2R53  \n')
	out.write('\t@angle:25  2.81862	127.8  0.00000	0.000    # NG2R53   CG2R53   OG2D1 \n')
	out.write('  }\n\n')
	out.write('   write_once("Data Dihedral Coeffs") {\n')
	out.write('\t@dihedral:1   0.01474	3	180	 0   # CG2R53   CG3C52   CG3C52   CG3C52  \n')
	out.write('\t@dihedral:2   0.00000	3	0    1   # CG2R53   CG3C52   CG3C52   HRA2  \n')
	out.write('\t@dihedral:3   0.03903	1	0    1   # CG2R53   NG2R53   CG311    CG321   \n')
	out.write('\t@dihedral:4   0.03903	1	0    1   # CG2R53   NG2R53   CG311    CG331   \n')
	out.write('\t@dihedral:5   0.00000	1	0    1   # CG2R53   NG2R53   CG311    HGA1  \n')
	out.write('\t@dihedral:6   0.10017	3	180	 0   # CG2R53   NG2R53   CG3C52   CG3C52  \n')
	out.write('\t@dihedral:7   0.00000	3	0    1   # CG2R53   NG2R53   CG3C52   HRA2  \n')
	out.write('\t@dihedral:8   0.00867	3	0    1   # CG311    CG321    CG311    CG321   \n')
	out.write('\t@dihedral:9   0.00867	3	0    1   # CG311    CG321    CG311    CG331   \n')
	out.write('\t@dihedral:10  0.00846	3	0    1   # CG311    CG321    CG311    HGA1  \n')
	out.write('\t@dihedral:11  0.00867	3	0    1   # CG311    CG321    CG311    NG2R53  \n')
	out.write('\t@dihedral:12  0.04857	2	180	 1   # CG311    NG2R53   CG2R53   CG3C52  \n')
	out.write('\t@dihedral:13  0.10841	2	180	 1   # CG311    NG2R53   CG2R53   OG2D1   \n')
	out.write('\t@dihedral:14  0.02905	3	0    1   # CG311    NG2R53   CG3C52   CG3C52  \n')
	out.write('\t@dihedral:15  0.00000	3	0    1   # CG311    NG2R53   CG3C52   HRA2  \n')
	out.write('\t@dihedral:16  0.00867	3	0    1   # CG321    CG311    CG321    HGA2  \n')
	out.write('\t@dihedral:17  0.00867	3	0    1   # CG321    CG311    CG331    HGA3  \n')
	out.write('\t@dihedral:18  0.00000	1	0    1   # CG321    CG311    NG2R53   CG3C52  \n')
	out.write('\t@dihedral:19  0.00867	3	0    1   # CG331    CG311    CG321    HGA2  \n')
	out.write('\t@dihedral:20  0.00000	1	0    1   # CG331    CG311    NG2R53   CG3C52  \n')
	out.write('\t@dihedral:21  0.01735	2	180	 0   # CG3C52   CG2R53   NG2R53   CG3C52  \n')
	out.write('\t@dihedral:22  0.04553	3	180	 0   # CG3C52   CG3C52   CG2R53   NG2R53  \n')
	out.write('\t@dihedral:23  0.00347	3	0    1   # CG3C52   CG3C52   CG2R53   OG2D1   \n')
	out.write('\t@dihedral:24  0.00824	3	0    1   # CG3C52   CG3C52   CG3C52   HRA2  \n')
	out.write('\t@dihedral:25  0.09236	3	0    0   # CG3C52   CG3C52   CG3C52   NG2R53  \n')
	out.write('\t@dihedral:26  0.11231	2	180	 1   # CG3C52   NG2R53   CG2R53   OG2D1   \n')
	out.write('\t@dihedral:27  0.00000	1	0    1   # CG3C52   NG2R53   CG311    HGA1  \n')
	out.write('\t@dihedral:28  0.00846	3	0    1   # HGA1 CG311    CG321    HGA2  \n')
	out.write('\t@dihedral:29  0.00846	3	0    1   # HGA1 CG311    CG331    HGA3  \n')
	out.write('\t@dihedral:30  0.00867	3	0    1   # HGA2 CG321    CG311    NG2R53    \n')
	out.write('\t@dihedral:31  0.00867	3	0    1   # HGA3 CG331    CG311    NG2R53    \n')
	out.write('\t@dihedral:32  0.00000	3	180	 1   # HRA2 CG3C52   CG2R53   NG2R53    \n')
	out.write('\t@dihedral:33  0.00000	3	0  	 1   # HRA2 CG3C52   CG2R53   OG2D1     \n')
	out.write('\t@dihedral:34  0.00824	3	0  	 1   # HRA2 CG3C52   CG3C52   HRA2  \n')
	out.write('\t@dihedral:35  0.00000	3	180	 1   # HRA2 CG3C52   CG3C52   NG2R53    \n')
	out.write('  }\n\n')
	out.write('   write_once("Data Improper Coeffs") {\n')
	out.write('\t@improper:1   3.90271	0.00000	 # CG2R53   CG3C52   NG2R53   OG2D1\n')
	out.write('  }\n\n')
  
	out.write('   write("Data Atoms") {\n')
	for i in range(len(atoms)):
		out.write('\t$atom:' + atoms[i][0] + ' \t\t$mol:.  @atom:' + atoms[i][2] + ' \t\t' + atoms[i][3] + ' \t\t' + atoms[i][4] + ' \t\t' + atoms[i][5] + ' \t\t' + atoms[i][6] + ' \t\t# ' + atoms[i][8] + '\n')
	out.write('  }\n\n   write("Data Bonds") {\n')
	for i in range(len(bonds)):
		out.write('\t$bond:' + bonds[i][0] + ' \t@bond:' + bonds[i][1] + ' \t$atom:' + bonds[i][2] + ' \t$atom:' + bonds[i][3] + ' \t# ' + bonds[i][5] + ' \t' + bonds[i][6] + '\n')
	out.write('  }\n\n   write("Data Angles") {\n')
	for i in range(len(angles)):
		out.write('\t$angle:' + angles[i][0] + ' \t@angle:' + angles[i][1] + ' \t$atom:' + angles[i][2] + ' \t$atom:' + angles[i][3] + ' \t$atom:' + angles[i][4] + ' \t# ' + angles[i][6] + ' \t' + angles[i][7] + ' \t' + angles[i][8] + '\n')
	out.write('  }\n\n   write("Data Dihedrals") {\n')
	for i in range(len(dihedrals)):
		out.write('\t$dihedral:' + dihedrals[i][0] + ' \t@dihedral:' + dihedrals[i][1] + ' \t$atom:' + dihedrals[i][2] + ' \t$atom:' + dihedrals[i][3] + ' \t$atom:' + dihedrals[i][4] + ' \t$atom:' + dihedrals[i][5] + ' \t# ' + dihedrals[i][7] + ' \t' + dihedrals[i][8] + ' \t' + dihedrals[i][9] + ' \t' + dihedrals[i][10] + '\n')
	out.write('  }\n\n   write("Data Impropers") {\n')
	for i in range(len(impropers)):
		out.write('\t$improper:' + impropers[i][0] + ' \t@improper:' + impropers[i][1] + ' \t$atom:' + impropers[i][2] + ' \t$atom:' + impropers[i][3] + ' \t$atom:' + impropers[i][4] + ' \t$atom:' + impropers[i][5] + ' \t# ' + impropers[i][7] + ' \t' + impropers[i][8] + ' \t' + impropers[i][9] + ' \t' + impropers[i][10] + '\n')
	out.write('  }\n\n}\n')