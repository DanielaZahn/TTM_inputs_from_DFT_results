import scipy.interpolate
from scipy.special import expit
import scipy.constants as constants
import numpy as np

def calculate_chemical_potential(temperatures,e_dos,fermi_energy,\
                               precision=0.000000001,initial_step=0.5,region_around_fermi_energy=5,number_of_points = 5000):
    # This function calculates the chemical potential as a function of temperature.
    # INPUTS:
    # temperatures: array of shape (N,) corresponding to the temperatures for which the chemical potential is calculated.
    # e_dos: array of shape (M,2) which contains information on the electronic density of states. The first column is 
    # energy in eV and the second column is the corresponding elecronic DOS in states/eV. 
    # Fermi energy: number corresponding to the Fermi energy, in eV. 
    # OPTIONAL INPUTS:
    # (Usually, the default values for all optional inputs should work, so don't change them without a reason.)
    # precision: number corresponding to the desired precision of the chemical potential in eV. If this value is too large, 
    # the calculations of G and especially c_e become very noisy or fail. If it is very small, the calculation becomes slow.
    # initial_step: number corresponding to the initial step size employed for the bisection. This only affects the efficiency 
    # of the bisection. It should not be a very small number. 
    # region_around_fermi_energy: number, defines the region around the Fermi energy that is considered in the calculation, in eV. 
    # The default value of 5 (eV) is more than enough to capture everything relevant. 
    # number_of_points: number, defines the number interpolation points of the electronic density of states in the region around 
    # the Fermi energy. This number should not be too small. 
    # 
    # OUTPUT:
    # mu: array of shape (N,) corresponding to the chemical potential at the temperatures given in "temperatures", in eV.
    

    # constants
    K_B = constants.Boltzmann*constants.physical_constants['joule-electron volt relationship'][0]  # Boltzmann constant, in eV/K
    
    # interpolate e_dos with very fine spacing in relevant region
    e_dos_function = scipy.interpolate.interp1d(e_dos[:,0],e_dos[:,1], kind='linear')
    e_dos = np.empty([number_of_points,2])
    # energies (fine spacing):
    e_dos[:,0] = np.linspace(fermi_energy-region_around_fermi_energy,fermi_energy+region_around_fermi_energy,number_of_points)
    # density of states (fine spacing): 
    e_dos[:,1] = e_dos_function(e_dos[:,0])

    # get number of electrons at T=0K
    n0 = np.trapz(e_dos[:,1]*np.heaviside(-(e_dos[:,0]-fermi_energy),0.5),e_dos[:,0])
    
    # find chemical potential with bisection method 
    mu = np.empty(np.shape(temperatures))
    
    for i in range(0,len(temperatures)):
        if temperatures[i] == 0:
            mu[i] = fermi_energy
        else:
            fermi_dirac_distribution = expit(-(e_dos[:,0]-fermi_energy)/(K_B*temperatures[i]))
            number_of_electrons = np.trapz(e_dos[:,1]*fermi_dirac_distribution,e_dos[:,0])
            
            # starting parameters:
            step = initial_step # initial step size for the calculation (bisection method), in eV
            direction = np.empty([1])
            not_yet_found = True # boolean that turns false when the chemical potential is found (within desired precision)
            mu[i] = fermi_energy # initial guess for chemical potential

            # bisection 
            while not_yet_found:
                if number_of_electrons < n0:
                    if direction == -1:
                        step = step/2
                    mu[i] = mu[i]+step
                    direction = 1
                elif number_of_electrons > n0:
                    if direction == 1:
                        step = step/2
                    mu[i] = mu[i]-step
                    direction = -1
                # calculate particle number with new chemical potential
                fermi_dirac_distribution = expit(-(e_dos[:,0]-mu[i])/(K_B*temperatures[i]))
                number_of_electrons = np.trapz(e_dos[:,1]*fermi_dirac_distribution,e_dos[:,0])
                if step < precision:
                    not_yet_found = False
    return mu