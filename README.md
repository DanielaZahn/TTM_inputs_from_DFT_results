# Calculation of the electron-phonon coupling parameter and heat capacities from DFT results (for metals)

This code calculates input parameters for a two-temperature model (TTM) or similar energy flow models from density functional theory (DFT) results.

Note that the calculation works only for metals. 

More precisely, the electron-phonon coupling parameter (G_ep), the electron heat capacity and the phonon heat capacity are calculated. 
The electron-phonon coupling parameter is calculated as in Waldecker et al., Phys. Rev. X 6, 021003 (https://doi.org/10.1103/PhysRevX.6.021003).

The required DFT results are: electron density of states, phonon density of states, Fermi energy, and Eliashberg function. 

There are two options depending on whether the DFT calculation is spin resolved: For spin-resolved DFT calculations, use "main_spin_resolved.ipynb". Otherwise, use "main_not_spin_resolved.ipynb".


To calculate TTM input parameters for your material, follow these four simple steps:

1) Copy your DFT results to the folder "inputs".
2) Write a python script that reads the DFT results and save this script in the folder "load_inputs". The name of the script needs to be the material name plus "_spin_resolved.py" or "_not_spin_resolved.py". Follow the examples provided in "load_inputs" for iron, cobalt, and nickel. Note that all variable names have to be the same as in the examples.
3) Change the file "config.cfg" such that the material specified there is your material. In addition, you can change other settings like lattice temperature (for calculation of G_ep) and the temperature range for which you would like the TTM inputs. 
4) Run the jupyter notebook "main_spin_resolved.ipynb" or "main_not_spin_resolved.ipynb", depending on which case applies for your DFT calculation. The results will be saved in the folder "results".


Additional notes:

- As an alternative to the jupyter notebooks, also a python script version of the main scripts is included ("main_spin_resolved.py" and "main_not_spin_resolved.py"). 
- It is also possible to calculate only the heat capacities and not G_ep, e.g. if you want to fit G_ep in your TTM. In this case, the Eliashberg function is not required. Simply don't run the cell that calculates G_ep in the jupyter notebook  / comment out the part of the code that calculates G_ep in the python script version. 
