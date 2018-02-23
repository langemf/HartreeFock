# The point of this file is to help you check your implementation along the way
# The first half builds your system, so similarly to the sample_driver.py file, change this as necessary
# For now, only the molecules that I've run are testable

import time
import mol
import tools
import scf
import integrals
import numpy as np
import pickle

# Specific a basis set that you would like to use (eg. sto-3g)
bas = 'sto3g'
bas = '631g'

# Give the name of the file for your molecule (eg. h2o.dat)
mol_name = 'co'
molec = mol_name
mol_name = 'test_molecules/' + mol_name # Change this depending on where your coordinate file is

# Specify the overall charge of the system (eg. -1,0,1)
charge = 0 

# Build the molecule Object
mole = mol.Molecule(bas, mol_name, charge)

# Create the SCF object for the chosen molecule
hf = scf.SCF(mole)

# If you're still building, it might be useful to build one matrix at a time
#For example, the overlap matrix
S = hf.create_S()

# Finally, run the scf procedure
hf.kernel()

# Now you can call whatever matrices or energies you would like from the SCF object
print("Total energy: {}".format(hf.e_tot))


# I would recommend not changing anything down here. This portion loads the assumed "correct" HF object
# and checks some important metrics with your implementation

with open('checks/' + molec + '_' + bas + '.pkl', 'rb') as f:
	hf_check = pickle.load(f)

print('Checking Overlap Matrix...')
print(np.allclose(hf.S, hf_check.S))
print('Checking Kinetic Energy Matrix...')
print(np.allclose(hf.T, hf_check.T))
print('Checking Enuc...')
print(np.allclose(hf.Vne, hf_check.Vne))
print('Checking eris...')
print(np.allclose(hf.eris, hf_check.eris))
print('Checking Density Matrix...')
print(np.allclose(hf.P, hf_check.P))
print('Checking Electronic Energy...')
print(np.allclose(hf.e_e, hf_check.e_e))
print('Checking Total Energy...')
print(np.allclose(hf.e_tot, hf_check.e_tot))


