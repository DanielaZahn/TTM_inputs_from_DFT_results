#!/usr/bin/env python
# coding: utf-8

# # Calculation of G_ep and heat capacities from DFT results (not spin resolved)
# 
# This notebook calculates the electron-phonon coupling parameter $G_\mathrm{ep}$ and the electron and phonon heat capacity from DFT results (electronic density of states, phonon density of states, Fermi level, Eliashberg function). 
# 
# Note that this is the version for a non-spin-resolved DFT calculation. Use the notebook main_spin_resolved for a spin-resolved DFT calculation (e.g. for ferromagnetic materials).
# 

# ## Load required modules, settings from config file, DFT results and other material parameters
# 
# In order to perform the calculation for another material, create a .py file named material+'_spin_resolved' in the folder 'load_inputs' that loads all necessary material-specific data (see examples). Furthermore, change the material entry in the file confic.cfg.


import numpy as np
get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib.pyplot as plt


# import functions saved separately
from functions.chemical_potential import calculate_chemical_potential
from functions.electron_phonon_coupling import calculate_electron_phonon_coupling
from functions.electron_heat_capacity import calculate_electron_heat_capacity
from functions.phonon_heat_capacity import calculate_phonon_heat_capacity

# import settings from config file
import configparser
config = configparser.ConfigParser()
config.read('config.cfg')
# material (string), to load the corresponding material-specific data and to save the results in the right text file:
material=config['GENERAL']['Material']
# lattice temperature in K, for the calculation of the electron-phonon coupling parameter G:
lattice_temperature=float(config['GENERAL']['Lattice_temperature'])
# desired temperature range in which all quantities should be calculated, in K:
temperatures = np.linspace(int(config['TEMPERATURE']['Temperature_min']),\
                           int(config['TEMPERATURE']['Temperature_max']),\
                           int(config['TEMPERATURE']['Temperature_points']))


# run the script that imports all required material-specific data (DFT calculation results and unit cell volume)
input_script_name='load_inputs/'+material+'_not_spin_resolved'
get_ipython().run_line_magic('run', '$input_script_name')


# ## Chemical potential
# The chemical potential is required for calculating the electron-phonon coupling and the electronic heat capacity.

mu = calculate_chemical_potential(temperatures,e_dos,fermi_energy)


# ## Electron-phonon coupling parameter ($G_\mathrm{ep}$)
# 
# The electron-phonon coupling parameter is calculated as in Waldecker et al., Phys. Rev. X 6, 021003 (https://doi.org/10.1103/PhysRevX.6.021003), Equation 9, with the (small) difference that changes of the chemical potential with electron temperature are considered here. 

g = calculate_electron_phonon_coupling(temperatures,lattice_temperature,mu,e_dos,fermi_energy,eliashberg)

# convert from W/(unit cell K) to W/(m^3 K)
g = g/unit_cell_volume

# save the result as text file
np.savetxt('results/'+material+'_notSpinResolved_electronPhononCoupling.txt',np.transpose([temperatures,g]),fmt='%07.2f %e',    header='material: '+material+' \n'
           'electron-phonon coupling parameter G_ep '
           '(DFT calculation not spin resolved) \n'
           'temperature [K], g_ep [W/(m^3K)] ')


# ## Electron heat capacity

electron_heat_capacity = calculate_electron_heat_capacity(temperatures,mu,e_dos,fermi_energy)
# Due to the numerical differentiation, the temperature sampling changes. Therefore, the output has two colums 
# with the first column corresponding to the new temperature points.

# convert from J/(unit cell K) to J/(m^3 K)
electron_heat_capacity[:,1] = electron_heat_capacity[:,1]/unit_cell_volume

# save the result as text file
np.savetxt('results/'+material+'_notSpinResolved_electronHeatCapacity.txt',electron_heat_capacity,fmt='%07.2f %e',    header='material: '+material+' \n'
           'electron heat capacity (DFT calculation not spin resolved) \n'
           'temperature [K], heat capacity [J/(m^3K)] ')


# ## Phonon heat capacity

phonon_heat_capacity = calculate_phonon_heat_capacity(temperatures,v_dos)
# Due to the numerical differentiation, the temperature sampling changes. Therefore, the output has two colums 
# with the first column corresponding to the new temperature points.

# convert from J/(unit cell K) to J/(m^3 K)
phonon_heat_capacity[:,1] = phonon_heat_capacity[:,1]/unit_cell_volume

# save the result as text file
np.savetxt('results/'+material+'_notSpinResolved_phononHeatCapacity.txt',phonon_heat_capacity,fmt='%07.2f %e',    header='material: '+material+' \n'
           'phonon heat capacity (DFT calculation not spin resolved) \n'
           'temperature [K], heat capacity [J/(m^3K)] ')

print('Calculation successful. See folder "results".')