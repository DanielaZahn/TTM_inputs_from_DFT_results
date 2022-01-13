import scipy.constants as constants
import numpy as np
import re

# required constants
HARTREE_TO_EV = constants.physical_constants['joule-electron volt relationship'][0]\
      /constants.physical_constants['joule-hartree relationship'][0] # conversion factor from Hartree to eV
AVOGADROS_NUMBER = constants.Avogadro # particles/mol

# material-specific data

# read electronic DOS (preferably per unit cell, see unit cell volume)
# units here: Hartree, states per Hartree per unit cell (is converted to eV below)
e_dos_majo = np.loadtxt('inputs/Fe_spinResolved_eDOS_majority.txt')
e_dos_mino = np.loadtxt('inputs/Fe_spinResolved_eDOS_minority.txt')

# unit cell volume (or, if e_dos and v_dos are not given per unit cell, corresponding other volume)
# Here, the unit cell volume is calculated from the molar volume. (this works for materials with 
# one atom per primitive unit cell)
molar_volume = 7.09e-6 # m^3/mol
unit_cell_volume = molar_volume/AVOGADROS_NUMBER # m^3 per unit cell
# IMPORTANT: The volume of the variable "unitCellVolume" has to match the units of the densities of states.
# Otherwise the heat capacities and G_ep will be WRONG!! (by a factor)
# For example, here e_dos and v_dos are in units of states per eV PER UNIT CELL and the corresponding volume
# is the unit cell volume.


# read Fermi energy (which is the same for both spin types of course)
file = open('inputs/Fe_spinResolved_eDOS_majority.txt')
alltext = file.read()
file.close()
# find the line of the text file in which the Fermi energy is written
index1 = alltext.find('Fermi energy')
index2 = alltext[index1:].find('\n')
# find the number in this line (which is the Fermi energy)
fermi_energy = float(np.squeeze(re.findall('[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?',\
                                        alltext[index1:index1+index2])))

# convert e_dos and fermi_energy from Hartree to eV
e_dos_majo[:,0] = e_dos_majo[:,0]*HARTREE_TO_EV # energy needs to be in eV
e_dos_mino[:,0] = e_dos_mino[:,0]*HARTREE_TO_EV # energy needs to be in eV
fermi_energy=fermi_energy*HARTREE_TO_EV 
e_dos_majo[:,1] = e_dos_majo[:,1]/HARTREE_TO_EV # DOS needs to be in states per eV
e_dos_mino[:,1] = e_dos_mino[:,1]/HARTREE_TO_EV # DOS needs to be in states per eV

# load Eliashberg function
eliashberg_majo = np.loadtxt('inputs/Fe_spinResolved_EliashbergFunction_majority.txt') 
eliashberg_mino = np.loadtxt('inputs/Fe_spinResolved_EliashbergFunction_minority.txt') 
# convert energy from Hartree to eV
eliashberg_majo[:,0] = eliashberg_majo[:,0]*HARTREE_TO_EV # energy needs to be in eV
eliashberg_mino[:,0] = eliashberg_mino[:,0]*HARTREE_TO_EV # energy needs to be in eV
# the second column (the Eliashberg function) has no units and therefore doesn't need to be converted

# load phonon density of states
v_dos=np.loadtxt('inputs/Fe_spinResolved_vDOS.txt')
# convert energy from Hartree to eV
v_dos[:,0] = v_dos[:,0]*HARTREE_TO_EV # energy needs to be in eV
v_dos[:,1] = v_dos[:,1]/HARTREE_TO_EV # DOS needs to be in states per eV
v_dos = v_dos[:,0:2]
# optional double-check: integrating the phonon DOS has to yield 3 times the atoms per unit cell 
# (here: 1 atom per unit cell / integral has to be 3)
#print(np.trapz(v_dos[:,1],v_dos[:,0]))


print('Material-specific data for iron has been loaded.')

