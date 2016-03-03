import os
print os.getcwd()

with open("out.colvars.traj") as file:
	line = file.readline()
	print line
	line = file.readline()
	with open("timeseries.txt","w") as out:
		while line:
			rows = line.split()
			if rows[0] != '#':
				out.write(rows[0] + ' ' + rows[1] + '\n')
			line = file.readline()
