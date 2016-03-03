import os

for i in range(1,102):
	if not(os.path.isfile(str(i) + '_window/out.colvars.traj')):
		print(str(i) + ' : not exist +')
		continue
	if os.stat(str(i) + '_window/out.colvars.traj').st_size < 1641804:
		print(str(i) + ' : unfinished *' + str(os.stat(str(i) + '_window/out.colvars.traj').st_size))