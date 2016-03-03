import os
print os.getcwd()

import shutil, errno

def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise
		
os.makedirs("ui")
		
copyanything("../../7ns/k_0.7/ui/data","ui/data")

for i in range(1,102):
    with open(str(i) + "_window/out.colvars.traj") as file:
        line = file.readline()
        print i
	line = file.readline()
	line = file.readline()
        with open("ui/data/" + str(i),"a") as out:
            while line:
                rows = line.split()
		if rows[0] != '#':
			out.write(rows[1] + ' 16.14 ' + rows[2] + '\n')
		line = file.readline()
