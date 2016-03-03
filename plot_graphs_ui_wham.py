import numpy as np
import matplotlib.pyplot as plt
import os

system_name = 'PVP5 Ag(111) 1_IC 3ns k_1.7'

if not os.path.exists('graph'):
    os.makedirs('graph')

pmf_ui = np.genfromtxt('fe_ui.xy', comments='#')
pmf_wham = np.genfromtxt('fe_wham.xy', comments='#')
hist = np.genfromtxt('global_histogram.xy', comments='#')

pmf_ui[:,1] = pmf_ui[:,1] - max(pmf_ui[:,1])
pmf_wham[:,1] = pmf_wham[:,1] - max(pmf_wham[:,1])

plt.plot(pmf_ui[:,0], pmf_ui[:,1], label='UI')
plt.plot(pmf_wham[:,0], pmf_wham[:,1], label='WHAM')
plt.xlabel(r'z ($\AA$)')
plt.ylabel('PMF (kcal/mol)')
plt.legend()
plt.title(system_name)
plt.savefig('graph/pmf.png', bbox_inches='tight')
plt.show()
plt.close()

plt.plot(hist[:,0], hist[:,1])
plt.xlabel(r'z ($\AA$)')
plt.ylabel('Histogram Count')
plt.title(system_name)
plt.savefig('graph/histogram.png', bbox_inches='tight')
plt.show()
plt.close()