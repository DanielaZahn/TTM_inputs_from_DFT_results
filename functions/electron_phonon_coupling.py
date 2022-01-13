import scipy.interpolate
from scipy.special import expit
import scipy.constants as constants
import numpy as np


def calculate_electron_phonon_coupling(temperatures,lattice_temperature,mu,e_dos,fermi_energy,eliashberg,\
                                       region_around_fermi_energy=5,number_of_points = 5000):
    # This function calculates the electron-phonon coupling parameter G_ep as a function of electron temperature, based on 
    # Waldecker et al., Phys. Rev. X 6, 021003 (2016), Eq. 9
    # INPUTS:
    # IMPORTANT: All inputs need to be in the right units! (see information below)
    # temperatures: array of shape (N,) corresponding to the temperatures for which G_ep is calculated, in K.
    # lattice_temperature: number, lattice temperature for which G_ep is calculated, in K.
    # mu: array of shape (N,) corresponding to the chemical potential in eV at the temperatures given in "temperatures".
    # mu can be calculated using the function "calculateChemicalPotential".
    # Note that a high precision of the chemical potential is necessary. The default precision of "calculateChemicalPotential" 
    # should be sufficient.
    # e_dos: array of shape (M,2) which contains information on the electronic density of states. The first column is 
    # energy in eV and the second column is the corresponding electronic density of states in states/eV. 
    # Fermi energy: number corresponding to the Fermi energy, in eV. 
    # eliashberg: array of shape (L,2) which contains information on the Eliashberg function. The first column is energy in eV
    # and the second column is the corresponding Eliashberg function (no units). 
    # OPTIONAL INPUTS:
    # (Usually, the default values for all optional inputs should work, so don't change them without a reason.)
    # region_around_fermi_energy: number, defines the region around the Fermi energy that is considered in the calculation, in eV. 
    # The default value of 5 (eV) is more than enough to capture everything relevant. 
    # number_of_points: number, defines the number interpolation points of the electronic density of states in the region around 
    # the Fermi energy. This number should not be too small. 
    # 
    # OUTPUT:
    # g: array of shape (N,) corresponding to the electron-phonon coupling parameter G_ep at the electron temperatures 
    # given in "temperatures". The units depend on the units of the electronic density of states: If it is in states/eV/unit cell,
    # g is in J/unit cell as well. In other words, the volume to which e_dos and g refer is the same (here: one unit cell). 
    
    
    
    # constants
    K_B = constants.Boltzmann*constants.physical_constants['joule-electron volt relationship'][0]  # Boltzmann constant, in eV/K
    ELEMENTARY_CHARGE = constants.e # elementary charge in C
    HBAR = constants.hbar*constants.physical_constants['joule-electron volt relationship'][0] # reduced Planck constant, in eVs

    first_integral = np.empty(np.shape(temperatures))
    second_integral = np.empty(np.shape(temperatures))
    z = np.empty(np.shape(temperatures))
    g = np.empty(np.shape(temperatures))
    
    # interpolate e_dos with very fine spacing in relevant region
    e_dos_function = scipy.interpolate.interp1d(e_dos[:,0],e_dos[:,1], kind='linear')
    e_dos = np.empty([number_of_points,2])
    # energies (fine spacing):
    e_dos[:,0] = np.linspace(fermi_energy-region_around_fermi_energy,fermi_energy+region_around_fermi_energy,number_of_points)
    # density of states (fine spacing): 
    e_dos[:,1] = e_dos_function(e_dos[:,0])

    # ensure that eliashberg has only positive energy entries (zero causes numerical problems and the Eliashberg function
    # should anyway be zero at zero energy)
    eliashberg = eliashberg[eliashberg[:,0]>0,:]

    # Bose-Einstein distribution (temperature: lattice temperature, at energies eliashberg[:,0])
    bose_einstein_distribution_lattice_temperature = 1./(np.exp(eliashberg[:,0]/(K_B*lattice_temperature))-1)

    for i in range(0,len(temperatures)):
        # get DOS at the chemical potential
        e_dos_at_chemical_potential =e_dos_function(mu[i])
        # prefactor of Waldecker et al., Phys. Rev. X 6, 021003 (2016), Eq. 9
        # Note that here N_c is not in the prefactor, because here we calculate g per unit cell (assuming that the e_dos
        # is per unit cell). Furthermore, there is an additional HBAR here because here we integrate over E, not omega.
        prefactor = -2*np.pi/e_dos_at_chemical_potential/HBAR 
        prefactor = prefactor*ELEMENTARY_CHARGE # to make sure output is in J/(unit cell K) instead of eV/(unit cell K)
        # Fermi-Dirac distribution (electron temperature), at energies e_dos[:,0]
        if temperatures[i] == 0:
            fermi_dirac_distribution = np.heaviside(-(e_dos[:,0]-mu[i]),0.5)
        else:
            fermi_dirac_distribution = expit(-(e_dos[:,0]-mu[i])/(K_B*temperatures[i]))
        # Bose-Einstein distribution (electron temperature), at energies eliashberg[:,0]
        if temperatures[i] == 0:
            bose_einstein_distribution_electron_temperature = np.zeros(np.shape(eliashberg[:,0]))
            # this is ok because it was ensured above that all energies of eliashberg are >0
        else:
            bose_einstein_distribution_electron_temperature = 1./(np.exp(eliashberg[:,0]/(K_B*temperatures[i]))-1)
        # calculate derivative of Fermi-Dirac distribution with respect to E
        derivative_of_fermi_dirac_distribution = np.empty([np.shape(e_dos)[0]-1,np.shape(e_dos)[1]])
        derivative_of_fermi_dirac_distribution[:,0] = (e_dos[1:,0]+e_dos[0:-1,0])/2 # energies
        derivative_of_fermi_dirac_distribution[:,1] = np.diff(fermi_dirac_distribution)/(e_dos[1,0]-e_dos[0,0]) # derivative
        # get e_dos at these new energies (they changed because we calculated a numerical derivative with np.diff)
        e_dos_at_new_energies = np.empty([np.shape(e_dos)[0]-1,np.shape(e_dos)[1]])
        e_dos_at_new_energies[:,0] = derivative_of_fermi_dirac_distribution[:,0]
        e_dos_at_new_energies[:,1] = e_dos_function(e_dos_at_new_energies[:,0])
        # first integral of Waldecker et al., Phys. Rev. X 6, 021003 (2016), Eq. 9
        # the integration is over E here, not over omega as in Eq. 9, which adds 1/HBAR to the prefactor
        first_integral[i] = np.trapz(eliashberg[:,0]**2*eliashberg[:,1]*\
                                  (bose_einstein_distribution_electron_temperature-\
                                   bose_einstein_distribution_lattice_temperature),eliashberg[:,0])
        # second integral of Waldecker et al., Phys. Rev. X 6, 021003 (2016), Eq. 9 (most of the electron-temperature 
        # dependence of g is due to this part)
        second_integral[i] = np.trapz(derivative_of_fermi_dirac_distribution[:,1]*e_dos_at_new_energies[:,1]**2,\
                                              derivative_of_fermi_dirac_distribution[:,0])
        z[i] = prefactor*first_integral[i]*second_integral[i] #in W/unit cell (assuming that the e_dos is per unit cell)
        g[i] = z[i]/(temperatures[i]-lattice_temperature) # in W/(unit cell K) (assuming that the e_dos is per unit cell)
    return g
