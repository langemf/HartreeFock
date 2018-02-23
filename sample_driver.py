# This file is an example driver file that you can modify to whatever your needs are
# This is a minimum working example, but depending on what you implement or need, you can add a lot

import time
import mol
import tools
import scf
import integrals
import numpy as np
import pickle

# Specific a basis set that you would like to use (eg. sto-3g)
bas = 'sto3g'

# Give the name of the file for your molecule (eg. h2o.dat)
mol_name = 'h2'
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


# Uncomment this if you'd like to save your HF Object with pickle

#with open('checks/' + molec + '_' + bas + '.pkl', 'wb') as f:
#	pickle.dump(hf, f, pickle.HIGHEST_PROTOCOL)



