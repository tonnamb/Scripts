
# coding: utf-8

# In[1]:

from math import sqrt

fname1 = "atac_10mer_HRA2.data"
fname2 = "EG.data"
fname3 = "MB3.data"
fname_out = "out-MB3.pair"


# In[2]:

class Data:
    
    def __init__(self,fname):
        
        self.n_Atoms = 0
        self.n_Bonds = 0
        self.n_Anglesa = 0
        self.n_Dihedrals = 0
        self.n_Impropers = 0
        self.n_AtomTypes = 0
        self.n_BondTypes = 0
        self.n_AngleTypes = 0
        self.n_DihedralTypes = 0
        self.n_ImproperTypes = 0
        self.xhi = 0
        self.yhi = 0
        self.zhi = 0
        
        self.Masses = []
        self.Pair_Coeffs = []
        self.Bond_Coeffs = []
        self.Angle_Coeffs = []
        self.Dihedral_Coeffs = []
        self.Improper_Coeffs = []
        self.Atoms = []
        self.Velocities = []
        self.Bonds = []
        self.Angles = []
        self.Dihedrals = []
        self.Impropers = []
        
        with open(fname) as f:
            i_line = 0;
            self.i_line_Velocities = 100000000000;
            self.i_line_Improper_Coeffs = 100000000000;
            self.i_line_Impropers = 100000000000;
            for line in f:

                i_line = i_line + 1

                if not (line == "\n"):
                    if (line.strip().split()[-1] == 'atoms'):
                        self.n_Atoms = int(line.strip().split()[0])
                    if (line.strip().split()[-1] == 'bonds'):
                        self.n_Bonds = int(line.strip().split()[0])
                    if (line.strip().split()[-1] == 'angles'):
                        self.n_Angles = int(line.strip().split()[0])
                    if (line.strip().split()[-1] == 'dihedrals'):
                        self.n_Dihedrals = int(line.strip().split()[0])
                    if (line.strip().split()[-1] == 'impropers'):
                        self.n_Impropers = int(line.strip().split()[0])
                    if (line.strip().split()[-1] == 'types'):
                        if (line.strip().split()[1] == 'atom'):
                            self.n_AtomTypes = int(line.strip().split()[0])
                        if (line.strip().split()[1] == 'bond'):
                            self.n_BondTypes = int(line.strip().split()[0])
                        if (line.strip().split()[1] == 'angle'):
                            self.n_AngleTypes = int(line.strip().split()[0])
                        if (line.strip().split()[1] == 'dihedral'):
                            self.n_DihedralTypes = int(line.strip().split()[0])
                        if (line.strip().split()[1] == 'improper'):
                            self.n_ImproperTypes = int(line.strip().split()[0])
                    if (line.strip().split()[-1] == 'xhi'):
                        self.xhi = line.strip().split()[1]
                    if (line.strip().split()[-1] == 'yhi'):
                        self.yhi = line.strip().split()[1]
                    if (line.strip().split()[-1] == 'zhi'):
                        self.zhi = line.strip().split()[1]

                if (line.strip() == "Masses"):
                    self.i_line_Masses = i_line
                if (line.strip() == "Pair Coeffs"):
                    self.i_line_Pair_Coeffs = i_line
                if (line.strip() == "Bond Coeffs"):
                    self.i_line_Bond_Coeffs = i_line
                if (line.strip() == "Angle Coeffs"):
                    self.i_line_Angle_Coeffs = i_line
                if (line.strip() == "Dihedral Coeffs"):
                    self.i_line_Dihedral_Coeffs = i_line
                if (line.strip() == "Improper Coeffs"):
                    self.i_line_Improper_Coeffs = i_line

                if (line.strip() == "Atoms"):
                    self.i_line_Atoms = i_line
                if (line.strip() == "Velocities"):
                    self.i_line_Velocities = i_line
                if (line.strip() == "Bonds"):
                    self.i_line_Bonds = i_line
                if (line.strip() == "Angles"):
                    self.i_line_Angles = i_line
                if (line.strip() == "Dihedrals"):
                    self.i_line_Dihedrals = i_line
                if (line.strip() == "Impropers"):
                    self.i_line_Impropers = i_line

        i_line = 0;
        with open(fname) as f:
            for line in f:

                i_line = i_line + 1

                if ((self.i_line_Masses+2) <= i_line <= (self.i_line_Masses+self.n_AtomTypes+1)):
                    self.Masses.append(line)
                if ((self.i_line_Pair_Coeffs+2) <= i_line <= (self.i_line_Pair_Coeffs+self.n_AtomTypes+1)):
                    self.Pair_Coeffs.append(line)
                if ((self.i_line_Bond_Coeffs+2) <= i_line <= (self.i_line_Bond_Coeffs+self.n_BondTypes+1)):
                    self.Bond_Coeffs.append(line)            
                if ((self.i_line_Angle_Coeffs+2) <= i_line <= (self.i_line_Angle_Coeffs+self.n_AngleTypes+1)):
                    self.Angle_Coeffs.append(line)        
                if ((self.i_line_Dihedral_Coeffs+2) <= i_line <= (self.i_line_Dihedral_Coeffs+self.n_DihedralTypes+1)):
                    self.Dihedral_Coeffs.append(line)            
                if ((self.i_line_Improper_Coeffs+2) <= i_line <= (self.i_line_Improper_Coeffs+self.n_ImproperTypes+1)):
                    self.Improper_Coeffs.append(line)

                if ((self.i_line_Atoms+2) <= i_line <= (self.i_line_Atoms+self.n_Atoms+1)):
                    self.Atoms.append(line)
                if ((self.i_line_Velocities+2) <= i_line <= (self.i_line_Velocities+self.n_Atoms+1)):
                    self.Velocities.append(line)                
                if ((self.i_line_Bonds+2) <= i_line <= (self.i_line_Bonds+self.n_Bonds+1)):
                    self.Bonds.append(line)
                if ((self.i_line_Angles+2) <= i_line <= (self.i_line_Angles+self.n_Angles+1)):
                    self.Angles.append(line)
                if ((self.i_line_Dihedrals+2) <= i_line <= (self.i_line_Dihedrals+self.n_Dihedrals+1)):
                    self.Dihedrals.append(line)
                if ((self.i_line_Impropers+2) <= i_line <= (self.i_line_Impropers+self.n_Impropers+1)):
                    self.Impropers.append(line)
    
    def __del__(self):
          class_name = self.__class__.__name__
          print class_name, "destroyed"


# In[3]:

data1 = Data(fname1)
data2 = Data(fname2)
data3 = Data(fname3)


# In[4]:

data1_mod = Data(fname1)
data2_mod = Data(fname2)
data3_mod = Data(fname3)


# In[5]:

data1_mod.Pair_Coeffs = []
for line in data1.Pair_Coeffs:
    PairType = int(line.split()[0])
    PairCoeff1 = float(line.split()[1])*0.0433634
    PairCoeff2 = float(line.split()[2])
    PairCoeff3 = float(line.split()[3])*0.0433634
    TheRest = ' '.join(line.split()[4:])
    data1_mod.Pair_Coeffs.append("%i %f %f %f %s \n" % (PairType,PairCoeff1,PairCoeff2,PairCoeff3,TheRest))


# In[6]:

data2_mod.Pair_Coeffs = []
for line in data2.Pair_Coeffs:
    PairType = int(line.split()[0])+data1.n_AtomTypes
    PairCoeff1 = float(line.split()[1])*0.0433634
    PairCoeff2 = float(line.split()[2])
    PairCoeff3 = float(line.split()[3])*0.0433634
    TheRest = ' '.join(line.split()[4:])
    data2_mod.Pair_Coeffs.append("%i %f %f %f %s \n" % (PairType,PairCoeff1,PairCoeff2,PairCoeff3,TheRest))


# In[7]:

data3_mod.Pair_Coeffs = []
for line in data3.Pair_Coeffs:
    PairType = int(line.split()[0])+17
    PairCoeff1 = float(line.split()[1])*0.0433634
    PairCoeff2 = float(line.split()[2])
    PairCoeff3 = float(line.split()[3])*0.0433634
    TheRest = ' '.join(line.split()[4:])
    data3_mod.Pair_Coeffs.append("%i %f %f %f %s \n" % (PairType,PairCoeff1,PairCoeff2,PairCoeff3,TheRest))


# In[8]:

Pair_Coeffs = data1_mod.Pair_Coeffs + data2_mod.Pair_Coeffs + data3_mod.Pair_Coeffs


# In[11]:

with open(fname_out,"w") as out:
    counter = -1
    for line1 in Pair_Coeffs:
        line1 = line1.split()
        counter = counter + 1
        for line2 in Pair_Coeffs[counter:]:
            line2 = line2.split()
            out.write("pair_coeff %s %s lj/charmm/coul/long %f %f %f %f %s %s \n" % (line1[0],line2[0],
                                                                               sqrt(float(line1[1])*float(line2[1])),
                                                                               0.5*(float(line1[2])+float(line2[2])),
                                                                               sqrt(float(line1[3])*float(line2[3])),
                                                                               0.5*(float(line1[4])+float(line2[4])),
                                                                               ' '.join(line1[5:]),
                                                                               ' '.join(line2[5:])))


# In[10]:



