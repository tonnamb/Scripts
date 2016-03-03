import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import butter, filtfilt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

cutoff = '120'
fname = 'diff_cutoff_{0}/diff_cutoff_{0}.txt'.format(cutoff)
system = 'PVP20 Ag(111)'

diff = np.loadtxt(fname)

meanrc = np.loadtxt('meanrc.txt')

if not os.path.exists('diff_cutoff_{0}_smooth'.format(cutoff)):
    os.makedirs('diff_cutoff_{0}_smooth'.format(cutoff))

def butter_lowpass(normal_cutoff, order=5):
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filtfilt(data, normal_cutoff, order=5):
    b, a = butter_lowpass(normal_cutoff, order=order)
    y = filtfilt(b, a, data)
    return y

# interpolate + smooth
xx = np.linspace(meanrc.min(), meanrc.max(), 1000)
itp = interp1d(meanrc, diff, kind='linear')
window_size, poly_order = 101, 1
yy_sg = savgol_filter(itp(xx), window_size, poly_order)

# for i in range(5):
	# yy_sg = np.append(yy_sg, yy_sg[-1])
	# xx = np.append(xx, xx[-1]+0.2)

normal_cutoff = 0.013
lowpass_order = 5
diff_smooth = butter_lowpass_filtfilt(yy_sg, normal_cutoff, lowpass_order)

print diff_smooth[-1]

plt.plot(meanrc, diff, ':', label='raw')
plt.plot(xx, yy_sg, '--', label='savgol')
plt.plot(xx, diff_smooth, '-', label='low pass')
plt.xlabel(r'z ($\AA$)')
plt.ylabel(r'Diffusion coefficient ($cm^2/s$)')
plt.title(system+' Cutoff = {0}'.format(cutoff))
plt.legend(frameon=False, loc=0)
plt.savefig('diff_cutoff_{0}_smooth/diff_profile.png'.format(cutoff), bbox_inches='tight')
plt.show()
plt.close()

np.savetxt('diff_cutoff_{0}_smooth/diff_cutoff_{0}.txt'.format(cutoff), diff)
np.savetxt('diff_cutoff_{0}_smooth/meanrc.txt'.format(cutoff), meanrc)
np.savetxt('diff_cutoff_{0}_smooth/diff_smooth.txt'.format(cutoff), diff_smooth)
np.savetxt('diff_cutoff_{0}_smooth/meanrc_smooth.txt'.format(cutoff), xx)

with open('diff_cutoff_{0}_smooth/smooth_param.txt'.format(cutoff), 'w') as f:
	f.write('savgol filter:\n')
	f.write('window_size = {0}\n'.format(window_size))
	f.write('poly_order = {0}\n\n'.format(poly_order))
	f.write('lowpass filter:\n')
	f.write('normal_cutoff = {0}\n'.format(normal_cutoff))
	f.write('lowpass_order = {0}\n'.format(lowpass_order))