import numpy as np
import scipy.constants as constants

def calculate_phonon_heat_capacity(temperatures,v_dos):
    # This function calculates the phonon heat capacity from the phonon density of states.
    # INPUTS:
    # IMPORTANT: All inputs need to be in the right units! (see information below)
    # temperatures: array of shape (N,) corresponding to the temperatures for the heat capacity calculation, in K.
    # v_dos: array of shape (M,2) which contains information on the phonon density of states. The first column is 
    # energy in eV and the second column is the corresponding phonon density of states in states/eV. 
    #
    # OUTPUT: 
    # phonon_heat_capacity: array of shape(N-1,2) which contains information on the phonon heat capacity. The first
    # column is temperature (the temperature sampling is slightly different due to the numerical differentiation).
    # The second column is heat capacity. The units depend on the units of the phonon density of states: 
    # If it is in states/eV/unit cell, the heat capacity is in J/unit cell. In other words, the volume to which
    # v_dos and phonon_heat_capacity refer is the same (here: one unit cell). 
    
    # constants
    K_B = constants.Boltzmann*constants.physical_constants['joule-electron volt relationship'][0]  # Boltzmann constant, in eV/K
    ELEMENTARY_CHARGE = constants.e # elementary charge in C
    
    # ensure that v_dos has only positive energy entries (zero energy causes numerical problems and the v_dos
    # should anyway be zero there.)
    v_dos = v_dos[v_dos[:,0]>0,:]
    
    # step 1: calculate electron energy as a function of temperature
    phonon_energy = np.empty(np.shape(temperatures))
    phonon_heat_capacity = np.empty([np.shape(temperatures)[0]-1,2])
    for i in range(0,len(temperatures)): 
        if temperatures[i] == 0:
            phonon_energy[i] = 0
        else:
            bose_einstein_distribution = 1./(np.exp(v_dos[:,0]/(K_B*temperatures[i]))-1)
            thermal_occupation_number = bose_einstein_distribution*v_dos[:,1] # occupation number=distribution function*DOS
            phonon_energy[i] = np.trapz(thermal_occupation_number*v_dos[:,0],v_dos[:,0]) # energy stored in phonon system       
    # convert phonon energy from eV to J:
    phonon_energy = phonon_energy*ELEMENTARY_CHARGE
    
    # step 2: calculate numerical derivative of phonon energy with respect to temperature to obtain heat capacity
    # temperatures:
    phonon_heat_capacity[:,0] = (temperatures[0:-1]+temperatures[1:])/2
    # heat capacity in J/(unit cell K) (assuming that the DOS is per unit cell):
    phonon_heat_capacity[:,1] = np.diff(phonon_energy)/(temperatures[1]-temperatures[0])
    
    return phonon_heat_capacity