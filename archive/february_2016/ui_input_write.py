import os
print os.getcwd()

if not os.path.exists("ui/data"):
	os.makedirs("ui/data")

for i in range(1,102):
    with open(str(i) + "_window/out.colvars.traj") as file:
        line = file.readline()
        print i
        line = file.readline()
        with open("ui/data/" + str(i),"w") as out:
            while line:
                rows = line.split()
		if rows[0] != '#':
			out.write(rows[1] + ' 23.06 ' + rows[2] + '\n')
		line = file.readline()
