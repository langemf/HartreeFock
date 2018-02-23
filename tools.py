# This file helps to parse the basis set files downloaded from the EMSL database
# Right now, it only parses basis sets in the NWChem format
# These procedures are called in the molecule initialization part of the driver file
import mol
from collections import defaultdict
import copy
#ang2bohr = 1.88973
ang2bohr = 1
# !! A note on the above. The coordinate files I worked with were in Bohr units, but most likely
# yours will be in Angstrom. Use the above variable to use whatever units you need 

l_dict = {'S': 0, 'P': 1, 'D': 2}
Z_dict = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18}
def create_atoms(mol_name, bas):
	atom_dict = create_mol(bas)
	coord_file = mol_name.lower() + '.dat'
	atom_list = []
	with open(coord_file) as fil:
		for lin in fil:
			c = lin.strip().split()	
			atom = c[0]
			xyz = mol.Point(c[1],c[2],c[3])
			Z = Z_dict[atom]
			new_atom = mol.Atom(atom, xyz * ang2bohr, Z)
			for orb in atom_dict[atom]:
				new_orb = copy.copy(orb)
				new_orb.xyz = xyz * ang2bohr
				new_atom.orbitals.append(new_orb)
			atom_list.append(new_atom)
	return atom_list
				
def create_mol(bas):
	atom_dict = defaultdict(list)
	file_to_open = 'basis/' + bas + '.dat'
	with open(file_to_open, 'r') as fil:
		name = ''
		l = ''
		pgo = []
		orb = ''
		for line in fil:
			if not line.strip().startswith('#') and not line.strip().startswith('END') and not len(line.strip()) == 0:
				if line[0].isalpha():
					if orb == 'SP' and name != '':
						atom_dict[name].append(mol.Orbital('SP-S', pgo[0], [0,0,0]))
						atom_dict[name].append(mol.Orbital('SP-P', pgo[1], [1,0,0]))
						atom_dict[name].append(mol.Orbital('SP-P', pgo[1], [0,1,0]))
						atom_dict[name].append(mol.Orbital('SP-P', pgo[1], [0,0,1]))
					elif orb != 'SP' and name != '':
						atom_dict[name].append(mol.Orbital(orb, pgo, l))
					name = line.strip().split()[0]
					orb = line.strip().split()[1]
					if orb == 'SP':
						pgo = [[],[]]
					else:					
#						l = np.ones(3)*l_dict.get(line.strip().split()[1])
						l = [0,0,0]
						pgo = []
				else:
					li = line.strip().split()
					if orb == 'SP':
						pgo[0].append(mol.PGO(float(li[0]),float(li[1])))
						pgo[1].append(mol.PGO(float(li[0]),float(li[2])))
					else:
						pgo.append(mol.PGO(float(li[0]),float(li[1])))
	return atom_dict



