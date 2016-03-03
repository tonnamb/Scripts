import mdtraj
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import cm
from matplotlib import mlab as ml

system_name = '2P PVP10 half half Ag(100) surface'

if not os.path.exists('heatmap'):
    os.makedirs('heatmap')

class traj:
    def __init__(self, read, pdb, agatom_array, nbin, skip):
        t = mdtraj.load(read, top=pdb, unit_set='real')
        self.unitcell_x = t.unitcell_lengths[0][0]
        self.unitcell_y = t.unitcell_lengths[0][1]
        self.unitcell_z = t.unitcell_lengths[0][2]
        self.n_frame = t.n_frames-skip
        topology = t.topology
        topagatom = agatom_array
        topagposition = t.xyz[:, topagatom, 2]
        self.meantopag = np.mean(topagposition)
        oxygenatom = [atom.index for atom in topology.atoms if (atom.name == '21')]
        self.opos = t.xyz[skip:, oxygenatom]
        self.opos_x = t.xyz[skip:, oxygenatom, 0].flatten()
        self.opos_y = t.xyz[skip:, oxygenatom, 1].flatten()
        self.opos_z = t.xyz[skip:, oxygenatom, 2].flatten() - self.meantopag
        print read
        print self.n_frame
        
        self.y, binEdges = np.histogram(np.ndarray.flatten(self.opos_z), nbin, density=False)
        self.dist_bin = binEdges[1]-binEdges[0] # Unit = nm
        binEdges = 10*binEdges # Convert to Angstroms
        self.bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        self.y = self.y/(self.unitcell_x*self.unitcell_y*self.dist_bin*self.n_frame) # Convert to density
        self.y = self.y*1e-3 # *10e-3 to convert from nm to \AA^3
    
    def xsec(self, range_min, range_max, xyz_index):
        # xyz_index = 0 for x
        # xyz_index = 1 for y
        # xyz_index = 2 for z
        pos = self.opos[:, :, xyz_index]
        self.idx_array = []
        self.opos_xsec = []
        for i in range(len(pos)):
            idx = NP.where( (pos[i] > range_min)*(pos[i] < range_max) )[0]
            self.idx_array.append(idx)
            self.opos_xsec.append(self.opos[i][idx])
        self.opos_xsec = NP.asarray(self.opos_xsec)

    def prep_heatmap(self, xaxis, yaxis):
        x = []
        y = []
        for i in range(len(self.opos_xsec)):
            x.append(self.opos_xsec[i][:, xaxis].flatten())
            y.append(self.opos_xsec[i][:, yaxis].flatten())
        x = NP.asarray([item for sublist in x for item in sublist])
        y = NP.asarray([item for sublist in y for item in sublist])
        return x, y

traj = traj('traj_npt.lammpstrj', 'traj_npt.pdb', range(512, 576)+range(704,768), 100, 800)
gridsize=40

np.savetxt('density_o_pvp.txt', np.stack((traj.bincenters, traj.y)))

plt.plot(traj.bincenters,traj.y, 'r-', lw=1.5)
plt.xlabel('z (Angstrom)', fontsize=14)
plt.ylabel(r'$\rho_O$ (#/$\AA^3$)', fontsize=14)
plt.title(system_name, fontsize=14)
plt.savefig('average_density_profile_o_pvp.png', bbox_inches='tight')
plt.show()
plt.close()

# y vs z
plt.subplot(111)
plt.hexbin(traj.opos_y, traj.opos_z, gridsize=gridsize, cmap=cm.jet)
plt.axis([traj.opos_y.min(), traj.opos_y.max(), traj.opos_z.min(), traj.opos_z.max()])
cb = plt.colorbar()
cb.set_label('Relative density of oxygen atoms', fontsize=14)
plt.title(system_name+': y vs z')
plt.xlabel(r'y ($\AA$)', fontsize=14)
plt.ylabel(r'z ($\AA$)', fontsize=14)
plt.ylim(0.0, 1.0)
plt.savefig('heatmap/yz.png', bbox_inches='tight')
plt.show()
plt.close()

# x vs y
plt.subplot(111)
plt.hexbin(traj.opos_x, traj.opos_y, gridsize=gridsize, cmap=cm.jet)
plt.axis([traj.opos_x.min(), traj.opos_x.max(), traj.opos_y.min(), traj.opos_y.max()])
cb = plt.colorbar()
cb.set_label('Relative density of oxygen atoms', fontsize=14)
plt.title(system_name+': x vs y')
plt.xlabel(r'x ($\AA$)', fontsize=14)
plt.ylabel(r'y ($\AA$)', fontsize=14)
plt.savefig('heatmap/xy.png', bbox_inches='tight')
plt.show()
plt.close()

# x vs z
plt.subplot(111)
plt.hexbin(traj.opos_x, traj.opos_z, gridsize=gridsize, cmap=cm.jet)
plt.axis([traj.opos_x.min(), traj.opos_x.max(), traj.opos_z.min(), traj.opos_z.max()])
cb = plt.colorbar()
cb.set_label('Relative density of oxygen atoms', fontsize=14)
plt.title(system_name+': x vs z')
plt.xlabel(r'x ($\AA$)', fontsize=14)
plt.ylabel(r'z ($\AA$)', fontsize=14)
plt.ylim(0.0, 1.0)
plt.savefig('heatmap/xz.png', bbox_inches='tight')
plt.show()
plt.close()