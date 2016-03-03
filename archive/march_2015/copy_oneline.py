import os
print os.getcwd()

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
			break