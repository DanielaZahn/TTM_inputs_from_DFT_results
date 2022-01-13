import scipy.interpolate
from scipy.special import expit
import scipy.constants as constants
import numpy as np

def calculate_electron_heat_capacity(temperatures,mu,e_dos,fermi_energy,\
                                     region_around_fermi_energy=5,number_of_points = 5000):
    # This function calculates the electron heat capacity.
    # INPUTS:
    # IMPORTANT: All inputs need to be in the right units! (see information below)
    # temperatures: array of shape (N,) corresponding to the temperatures for which the chemical potential is calculated, in K.
    # mu: array of shape (N,) corresponding to the chemical potential at the temperatures given in "temperatures", in eV.
    # mu can be calculated using the function "calculateChemicalPotential".
    # Note that a very high precision of the chemical potential is necessary. The default precision of 
    # calculateChemicalPotential" should be sufficient.
    # e_dos: array of shape (M,2) which contains information on the electronic density of states. The first column is 
    # energy in eV and the second column is the corresponding electronic density of states in states/eV. 
    # fermi_energy: number corresponding to the Fermi energy, in eV. 
    # OPTIONAL INPUTS:
    # (Usually, the default values for all optional inputs should work, so don't change them without a reason.)
    # region_around_fermi_energy: number, defines the region around the Fermi energy that is considered in the calculation, in eV. 
    # The default value of 5 (eV) is more than enough to capture everything relevant. 
    # number_of_points: number, defines the number interpolation points of the electronic density of states in the region
    # around the Fermi energy. This number should not be too small. 
    #
    # OUTPUT: 
    # electron_heat_capacity: array of shape(N-1,2) which contains information on the electron heat capacity. The first
    # column is temperature (the temperature sampling is slightly different due to the numerical differentiation).
    # The second column is heat capacity. The units depend on the units of the electronic density of states: If it
    # is in states/eV/unit cell, the heat capacity is in J/unit cell. In other words, the volume to which e_dos and 
    # electron_heat_capacity refer is the same (here: one unit cell). 
    
    # constants
    K_B = constants.Boltzmann*constants.physical_constants['joule-electron volt relationship'][0]  # Boltzmann constant, in eV/K
    ELEMENTARY_CHARGE = constants.e # elementary charge in C
    
    # interpolate e_dos with very fine spacing in relevant region
    e_dos_function = scipy.interpolate.interp1d(e_dos[:,0],e_dos[:,1], kind='linear')
    e_dos = np.empty([number_of_points,2])
    # energies (fine spacing):
    e_dos[:,0] = np.linspace(fermi_energy-region_around_fermi_energy,fermi_energy+region_around_fermi_energy,number_of_points)
    # density of states (fine spacing): 
    e_dos[:,1] = e_dos_function(e_dos[:,0])
    
    # step 1: calculate electron energy 
    # it doesn't matter that only a window around E_f is considered here, 
    # because what matters for the electronic heat capacity are changes in energy, 
    # which only happen around E_f(or the chemical potential, to be precise, which is
    # close to E_f)
    energy=np.empty(np.shape(temperatures))
    for i in range(0,len(temperatures)):
        if temperatures[i] == 0:
            Fermi_Dirac_distribution=np.heaviside(-(e_dos[:,0]-mu[i]),0.5)
        else:
            Fermi_Dirac_distribution = expit(-(e_dos[:,0]-mu[i])/(K_B*temperatures[i]))
        energy[i] = np.trapz(e_dos[:,1]*Fermi_Dirac_distribution*e_dos[:,0],e_dos[:,0])
    # convert energy from eV to J:
    energy = energy*ELEMENTARY_CHARGE
   
    # step 2: calculate the heat capacity by differentiating energy 
    electron_heat_capacity = np.empty([np.size(temperatures,0)-1,2])
    # temperature:
    electron_heat_capacity[:,0] = (temperatures[0:-1]+temperatures[1:])/2
    # heat capacity in J/unit cell/K (assuming that the DOS is per unit cell):
    electron_heat_capacity[:,1] = np.diff(energy)/(temperatures[1]-temperatures[0])
    return electron_heat_capacity