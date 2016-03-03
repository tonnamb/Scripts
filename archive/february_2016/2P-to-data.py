
# coding: utf-8

# In[10]:

fname_2P = "MB3.data"
fname_del = "del.data"
fname_out = "out-MB3.data"


# In[11]:

class Data:
    
    def __init__(self,fname):
        
        self.n_Atoms = 0
        self.n_Bonds = 0
        self.n_Angles = 0
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


# In[12]:

data_2P = Data(fname_2P)
data_del = Data(fname_del)


# In[13]:

mol_id_del = [];
for line in data_del.Atoms:
    mol_id_del.append(int(line.split()[1]))
max_mol_id_del = max(mol_id_del)


# In[14]:

data_2P_mod = Data(fname_2P)


# In[15]:

data_2P_mod.Masses = []
for line in data_2P.Masses:
    AtomType = int(line.split()[0]) + data_del.n_AtomTypes
    TheRest = ' '.join(line.split()[1:])
    data_2P_mod.Masses.append("%i %s \n" % (AtomType,TheRest))


# In[16]:

data_2P_mod.Bond_Coeffs = []
for line in data_2P.Bond_Coeffs:
    BondType = int(line.split()[0]) + data_del.n_BondTypes
    TheRest = ' '.join(line.split()[1:])
    data_2P_mod.Bond_Coeffs.append("%i %s \n" % (BondType,TheRest))


# In[17]:

data_2P_mod.Angle_Coeffs = []
for line in data_2P.Angle_Coeffs:
    AngleType = int(line.split()[0]) + data_del.n_AngleTypes
    TheRest = ' '.join(line.split()[1:])
    data_2P_mod.Angle_Coeffs.append("%i %s \n" % (AngleType,TheRest))


# In[18]:

data_2P_mod.Dihedral_Coeffs = []
for line in data_2P.Dihedral_Coeffs:
    DihedralType = int(line.split()[0]) + data_del.n_DihedralTypes
    TheRest = ' '.join(line.split()[1:])
    data_2P_mod.Dihedral_Coeffs.append("%i %s \n" % (DihedralType,TheRest))


# In[19]:

data_2P_mod.Improper_Coeffs = []
for line in data_2P.Improper_Coeffs:
    ImproperType = int(line.split()[0]) + data_del.n_ImproperTypes
    TheRest = ' '.join(line.split()[1:])
    data_2P_mod.Improper_Coeffs.append("%i %s \n" % (ImproperType,TheRest))


# In[20]:

data_2P_mod.Atoms = []
data_2P_mod.Velocities = []
for line in data_2P.Atoms:
    AtomID = int(line.split()[0]) + data_del.n_Atoms
    MolID = int(line.split()[1]) + max_mol_id_del
    AtomType = int(line.split()[2]) + data_del.n_AtomTypes
    Coord = ' '.join(line.split()[3:7])
    TheRest = ' '.join(line.split()[7:])
    data_2P_mod.Atoms.append("%i %i %i %s 0 0 0 %s \n" % (AtomID,MolID,AtomType,Coord,TheRest))
    data_2P_mod.Velocities.append("%i 0 0 0 \n" % AtomID)


# In[21]:

data_2P_mod.Bonds = []
for line in data_2P.Bonds:
    BondID = int(line.split()[0]) + data_del.n_Bonds
    BondType = int(line.split()[1]) + data_del.n_BondTypes
    AtomsID1 = int(line.split()[2]) + data_del.n_Atoms
    AtomsID2 = int(line.split()[3]) + data_del.n_Atoms
    TheRest = ' '.join(line.split()[4:])
    data_2P_mod.Bonds.append("%i %i %i %i %s \n" % (BondID,BondType,AtomsID1,AtomsID2,TheRest))


# In[22]:

data_2P_mod.Angles = []
for line in data_2P.Angles:
    AngleID = int(line.split()[0]) + data_del.n_Angles
    AngleType = int(line.split()[1]) + data_del.n_AngleTypes
    AtomsID1 = int(line.split()[2]) + data_del.n_Atoms
    AtomsID2 = int(line.split()[3]) + data_del.n_Atoms
    AtomsID3 = int(line.split()[4]) + data_del.n_Atoms
    TheRest = ' '.join(line.split()[5:])
    data_2P_mod.Angles.append("%i %i %i %i %i %s \n" % (AngleID,AngleType,AtomsID1,AtomsID2,AtomsID3,TheRest))


# In[23]:

data_2P_mod.Dihedrals = []
for line in data_2P.Dihedrals:
    DihedralID = int(line.split()[0]) + data_del.n_Dihedrals
    DihedralType = int(line.split()[1]) + data_del.n_DihedralTypes
    AtomsID1 = int(line.split()[2]) + data_del.n_Atoms
    AtomsID2 = int(line.split()[3]) + data_del.n_Atoms
    AtomsID3 = int(line.split()[4]) + data_del.n_Atoms
    AtomsID4 = int(line.split()[5]) + data_del.n_Atoms
    TheRest = ' '.join(line.split()[6:])
    data_2P_mod.Dihedrals.append("%i %i %i %i %i %i %s \n" % (DihedralID,DihedralType,AtomsID1,AtomsID2,AtomsID3,AtomsID4,TheRest))


# In[24]:

data_2P_mod.Impropers = []
for line in data_2P.Impropers:
    ImproperID = int(line.split()[0]) + data_del.n_Impropers
    ImproperType = int(line.split()[1]) + data_del.n_ImproperTypes
    AtomsID1 = int(line.split()[2]) + data_del.n_Atoms
    AtomsID2 = int(line.split()[3]) + data_del.n_Atoms
    AtomsID3 = int(line.split()[4]) + data_del.n_Atoms
    AtomsID4 = int(line.split()[5]) + data_del.n_Atoms
    TheRest = ' '.join(line.split()[6:])
    data_2P_mod.Impropers.append("%i %i %i %i %i %i %s \n" % (ImproperID,ImproperType,AtomsID1,AtomsID2,AtomsID3,AtomsID4,TheRest))


# In[25]:

data_2P_mod.n_Atoms = data_2P.n_Atoms + data_del.n_Atoms
data_2P_mod.n_Bonds = data_2P.n_Bonds + data_del.n_Bonds
data_2P_mod.n_Angles = data_2P.n_Angles + data_del.n_Angles
data_2P_mod.n_Dihedrals = data_2P.n_Dihedrals + data_del.n_Dihedrals
data_2P_mod.n_Impropers = data_2P.n_Impropers + data_del.n_Impropers
data_2P_mod.n_AtomTypes = data_2P.n_AtomTypes + data_del.n_AtomTypes
data_2P_mod.n_BondTypes = data_2P.n_BondTypes + data_del.n_BondTypes
data_2P_mod.n_AngleTypes = data_2P.n_AngleTypes + data_del.n_AngleTypes
data_2P_mod.n_DihedralTypes = data_2P.n_DihedralTypes + data_del.n_DihedralTypes
data_2P_mod.n_ImproperTypes = data_2P.n_ImproperTypes + data_del.n_ImproperTypes


# In[26]:

with open(fname_out,"w") as out:
    out.write("LAMMPS Description # generated from combining %s and %s with 2P-to-data.ipynb \n\n" % (fname_2P,fname_del))
    out.write("%i atoms\n" % data_2P_mod.n_Atoms)
    out.write("%i atom types\n" % data_2P_mod.n_AtomTypes)
    out.write("%i bonds\n" % data_2P_mod.n_Bonds)
    out.write("%i bond types\n" % data_2P_mod.n_BondTypes)
    out.write("%i angles\n" % data_2P_mod.n_Angles)
    out.write("%i angle types\n" % data_2P_mod.n_AngleTypes)
    out.write("%i dihedrals\n" % data_2P_mod.n_Dihedrals)
    out.write("%i dihedral types\n" % data_2P_mod.n_DihedralTypes)
    out.write("%i impropers\n" % data_2P_mod.n_Impropers)
    out.write("%i improper types\n\n" % data_2P_mod.n_ImproperTypes)
    out.write("0.0000000000000000e+00 %s xlo xhi\n" % data_del.xhi)
    out.write("0.0000000000000000e+00 %s ylo yhi\n" % data_del.yhi)
    out.write("0.0000000000000000e+00 %s zlo zhi\n\n" % data_del.zhi)
    out.write("Masses\n\n")
    for line in data_del.Masses:
        out.write(line)
    for line in data_2P_mod.Masses:
        out.write(line)
    out.write("\nBond Coeffs\n\n")
    for line in data_del.Bond_Coeffs:
        out.write(line)
    for line in data_2P_mod.Bond_Coeffs:
        out.write(line)
    out.write("\nAngle Coeffs\n\n")
    for line in data_del.Angle_Coeffs:
        out.write(line)
    for line in data_2P_mod.Angle_Coeffs:
        out.write(line)
    out.write("\nDihedral Coeffs\n\n")
    for line in data_del.Dihedral_Coeffs:
        out.write(line)
    for line in data_2P_mod.Dihedral_Coeffs:
        out.write(line)
    out.write("\nImproper Coeffs\n\n")
    for line in data_del.Improper_Coeffs:
        out.write(line)
    for line in data_2P_mod.Improper_Coeffs:
        out.write(line)
    out.write("\nAtoms\n\n")
    for line in data_del.Atoms:
        out.write(line)
    for line in data_2P_mod.Atoms:
        out.write(line)
    out.write("\nVelocities\n\n")
    for line in data_del.Velocities:
        out.write(line)
    for line in data_2P_mod.Velocities:
        out.write(line)
    out.write("\nBonds\n\n")
    for line in data_del.Bonds:
        out.write(line)
    for line in data_2P_mod.Bonds:
        out.write(line)
    out.write("\nAngles\n\n")
    for line in data_del.Angles:
        out.write(line)
    for line in data_2P_mod.Angles:
        out.write(line)
    out.write("\nDihedrals\n\n")
    for line in data_del.Dihedrals:
        out.write(line)
    for line in data_2P_mod.Dihedrals:
        out.write(line)
    out.write("\nImpropers\n\n")
    for line in data_del.Impropers:
        out.write(line)
    for line in data_2P_mod.Impropers:
        out.write(line)


# In[26]:



