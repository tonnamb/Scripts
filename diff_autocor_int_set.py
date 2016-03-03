import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import os

system_name = r'PVP5 Ag(111) k = 1.7 eV/$\AA^2$'
spec_folder = "data"
spec_nwin = 60
cutoff_array = [30, 60, 90, 120, 150, 200]

for i in range(5):
    for c in range(len(cutoff_array)):
        cutoff = cutoff_array[c]
        if not os.path.exists('set_{0}/diff_cutoff_{1}'.format(i+1, cutoff)):
            os.makedirs('set_{0}/diff_cutoff_{1}'.format(i+1, cutoff))

for c in range(len(cutoff_array)):
    cutoff = cutoff_array[c]
    if not os.path.exists('set_mean/diff_cutoff_{0}'.format(cutoff)):
        os.makedirs('set_mean/diff_cutoff_{0}'.format(cutoff))

def estimated_autocorrelation(x):
	# http://stackoverflow.com/questions/14297012/estimate-autocorrelation-using-python
    n = len(x)
    variance = np.var(x)
    x = x-np.mean(x)
    r = np.correlate(x, x, mode = 'full')[-n:]
    assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
    result = r/(variance*(np.arange(n, 0, -1)))
    return result
            
class Data:

    def __init__(self, folder, nwin, start, end):
        self.RC = [[] for i in range(nwin)]
        self.nwin = nwin
        self.Center = []
        self.MeanRC = []
        self.VarRC = []
        self.NStep = []

        for i in range(nwin):
            with open(folder + '/' + str(i+1)) as f:
                run_once = 0
                count = 0
                for line in f:
                    if (count >= start):
                        self.RC[i].append(float(line.strip().split()[0]))
                    count = count + 1
                    if run_once == 0:
                        self.Center.append(float(line.strip().split()[2]))
                        self.Spring = float(line.strip().split()[1])
                        run_once = 1
                    if (count > end):
                        break
            self.MeanRC.append(np.mean(self.RC[i]))
            self.VarRC.append(np.var(self.RC[i]))
            self.NStep.append(len(self.RC[i]))

    def calculate_diffusion(self, cutoff_array, set_index):
        
        np.savetxt('set_{0}/meanrc.txt'.format(set_index), np.asarray(self.MeanRC))
        
        self.calc_diff_array = [[] for i in range(len(cutoff_array))]

        for i in range(self.nwin):
            x = range(self.NStep[i])
            y = estimated_autocorrelation(self.RC[i])
            for c in range(len(cutoff_array)):
                cutoff = cutoff_array[c]
                xcut = np.asarray(x[:cutoff])*0.0015*100
                # Conversion from timestep to ps [=] 0.0015 picosecond timesteps * 100 colvarsTrajFrequency
                ycut = y[:cutoff]
                plt.plot(xcut, ycut)
                plt.xlabel('time (ps)')
                plt.ylabel('autocorrelation')
                plt.title('Window {0}'.format(i+1))
                plt.savefig('set_{0}/diff_cutoff_{1}/auto_{2}.png'.format(set_index, cutoff, i+1), bbox_inches='tight')
                plt.close()

                integrate = simps(ycut, dx=150*10**(-15))
                # dx = 1 timestep = 15*10**(-15) seconds [=] 1.5 femtosecond * 100 colvarsTrajFrequency

                calc_diff = self.VarRC[i]*(10**(-8))**2/integrate
                # calc_diff [=] cm^2/s
                # (10**(-8))**2 is the conversion from Angstrom^2 to cm^2

                print 'Set {0}, Cutoff {1}, Window {2}: D = {3} cm^2/s'.format(set_index, cutoff, i+1, calc_diff)
                
                self.calc_diff_array[c].append(calc_diff)

        self.calc_diff_array = np.asarray(self.calc_diff_array)
        
        for c in range(len(cutoff_array)):
            cutoff = cutoff_array[c]
            np.savetxt('set_{0}/diff_cutoff_{1}/diff_cutoff_{1}.txt'.format(set_index, cutoff), self.calc_diff_array[c])
            plt.plot(self.MeanRC, self.calc_diff_array[c], label='{0} ({1} ps)'.format(cutoff_array[c], cutoff_array[c]*0.15))
        plt.ylabel(r'Diffusion coefficient ($cm^2/s$)')
        plt.xlabel(r'z ($\AA$)')
        plt.title(system_name)
        plt.legend(frameon=False, loc=0)
        plt.savefig('set_{0}/plot_diff_cutoff.png'.format(set_index), bbox_inches='tight')
        plt.close()

# if start=1, will start reading row 2
# if end=1, will read row 1 and 2
data_1 = Data(spec_folder, spec_nwin, 0, 19999)
data_2 = Data(spec_folder, spec_nwin, 20000, 29999)
data_3 = Data(spec_folder, spec_nwin, 30000, 39999)

data_1.calculate_diffusion(cutoff_array, 1)
data_2.calculate_diffusion(cutoff_array, 2)
data_3.calculate_diffusion(cutoff_array, 3)

mean_diff_array = []
for c in range(len(cutoff_array)):
    mean_diff_array.append((data_1.calc_diff_array[c] + data_2.calc_diff_array[c] + data_3.calc_diff_array[c])/3.0)
    set_meanrc = (np.asarray(data_1.MeanRC) + np.asarray(data_2.MeanRC) + np.asarray(data_3.MeanRC))/3.0

np.savetxt('set_mean/meanrc.txt', set_meanrc)

for c in range(len(cutoff_array)):
    cutoff = cutoff_array[c]
    np.savetxt('set_mean/diff_cutoff_{0}/diff_cutoff_{0}.txt'.format(cutoff), mean_diff_array[c])
    plt.plot(set_meanrc, mean_diff_array[c], label='{0} ({1} ps)'.format(cutoff_array[c], cutoff_array[c]*0.15))
plt.ylabel(r'Diffusion coefficient ($cm^2/s$)')
plt.xlabel(r'z ($\AA$)')
plt.title(system_name)
plt.legend(frameon=False, loc=0)
plt.savefig('set_mean/plot_diff_cutoff.png', bbox_inches='tight')
plt.close()