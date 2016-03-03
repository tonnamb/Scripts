import numpy as np
import matplotlib.pyplot as plt
import os
# print os.getcwd()

if not os.path.exists('timeseries'):
    os.makedirs('timeseries')

z = []
step = []
with open("out.colvars.traj") as file:
    line = file.readline()
    # print line
    line = file.readline()
    with open("timeseries/timeseries.txt","w") as out:
        while line:
            rows = line.split()
            if rows[0] != '#':
                out.write(rows[0] + ' ' + rows[1] + '\n')
                step.append(int(rows[0]))
                z.append(float(rows[1]))
            line = file.readline()

plt.plot(step, z)
plt.xlabel('Timesteps')
plt.ylabel(r'z ($\AA$)')
plt.savefig('timeseries/plot_timeseries.png', bbox_inches='tight')
plt.show()

diff = []
print 'Count of non-monotonicity'
for i in range(10):
    datatest = z[i::10]
    diff_temp = sum(np.diff(datatest)>=0)
    diff.append(diff_temp)
    print 'From {0}: {1}'.format(step[i], diff_temp)

min_index = diff.index(min(diff))
print 'Choosing from: {0}'.format(step[min_index])

run_z = z[min_index::10][1:]
run_step = step[min_index::10][1:]
print run_step
bias_z = np.arange(29-0.45, 2-0.45, -0.45)
plt.scatter(run_step, run_z, c='r', label='Initial')
plt.scatter(run_step, bias_z, c='b', alpha=0.5, label='Bias')
plt.legend()
plt.xlabel('Timesteps')
plt.ylabel(r'z ($\AA$)')
plt.savefig('timeseries/compare_z.png', bbox_inches='tight')
plt.show()

print 'As start in write_1ns.py: {0}'.format(run_step[-1])
with open('timeseries/start_in_write1ns.txt', 'w') as f:
    f.write('As start in write_1ns.py: {0}'.format(run_step[-1]))