import os

os.makedirs("data")

with open("out.colvars.traj") as file:
    line = file.readline()
    line = file.readline()
    with open("data/1","w") as out:
        while line:
            rows = line.split()
	if rows[0] != '#':
		out.write(rows[1] + ' 16.14 ' + rows[2] + '\n')
	line = file.readline()