import numpy as np
import tools

# Molecule Object
# This object contains all atoms objects and the overall name and charge of the molecule of interest
# Additionally, all of the matrices involved in the SCF procedure are currently available in the Molecule object

class Molecule:

	def __init__(self, basis, name, charge):
		self.atoms = tools.create_atoms(name, basis)	
		self.name = name
		self.charge = charge 
		self.elec = self.elec_count()
	
	def __repr__(self):
		s = "Molecule(%s)" % (self.name)
		return s

	def __str__(self):
		s = "Molecule Object: %s" % (self.name)
		return s

	def nmo(self):
		count = 0
		for a in self.atoms:
			count += len(a.orbitals)
		return count

	def elec_count(self):
		count = 0
		for a in self.atoms:
			count += a.Z
		return count - self.charge


# Atom class which represents the individual atoms in a molecule object
# This includes the name, coordinates, and charge of the atom
# as well as the atomic orbitals associated with the atom
# __repr__ and __str__ are implemented to help in potential logging needs
class Atom:
	
	def __init__(self):
		self.name = None
		self.orbitals = None
		self.xyz = None
		self.Z = None

	def __init__(self, name, xyz, Z):
		self.name = name
		self.orbitals = []
		self.xyz = xyz
		self.Z = Z
	'''
	def __init__(self, name, orbitals, xyz):
		self.name = name
		self.orbitals = orbital
		self.xyz = xyz
	'''
	def __repr__(self):
		s = "Atom(%s, %s, %s)" % (self.name, self.xyz, self.Z)
		return s	

	def __str__(self):
		s = "Atom Object\nAtom: %s \nXYZ: %s\nZ: %s" %(self.name, self.xyz, self.Z)
		return s

	def distance(self, other):
		return np.sqrt((self.xyz.x-other.xyz.x)**2+(self.xyz.y-other.xyz.y)**2+(self.xyz.z-other.xyz.z)**2)

# Orbital object which is used for atomic orbitals to keep track of position and PGOs
# The Orbital object is a class variable in the Atom object
# This helps keep track of important properties such as angular momentum (l), the position (xyz)
# the type (typ) and the PGO's making up the atomic orbital (pgo)
# __repr__ and __str__ are implemented to help in potential logging needs
class Orbital:

	def __init__(self):
		self.typ = None
		self.pgo = None
		self.l = None
		self.xyz = None

	def __init__(self, typ, pgo, l):
		self.typ = typ
		self.pgo = pgo
		self.l = l
		self.xyz = None
	
	def __repr__(self):
		s = "Orbital(%s, %s, %s, %s)" % (self.typ, self.xyz, self.l, self.pgo)
		return s	

	def __str__(self):
		s = "Orbital Object\nType: %s \nXYZ: %s\nl: %s \npgo: %s" %(self.typ, self.xyz, self.l, self.pgo)
		return s
# The PGO object which is a class variable in the Orbital object
# The PGO object houses the necessary exponents (exp) and coefficients (coeff) for the PGOs that make up a
# given atomic orbital
# __repr__ and __str__ are implemented to help in potential logging needs
class PGO:

	def __init__(self, exp, coeff):
		self.exp = None
		self.coeff = None
	
	def __init__(self, exp, coeff):
		self.exp = exp
		self.coeff = coeff
	
	def __repr__(self):
		s = "PGO(%s, %s)" % (self.exp, self.coeff)
		return s

	def __string__(self):
		s = "PGO Object\nExponent: %s\n,Coefficient: %s" % (self.exp, self.coeff)
		return s

# The Point object, which is probably overkill, but was great practice for implementing a variety of special methods
# The Point object is what it sounds like, it keeps track of the different coordinates of position which is used for 
# all xyz variables in the objects above.
# Many special methods were implemented for this object, mainly for practice
class Point:
	
	def __init__(self):
		self.x = None
		self.y = None
		self.z = None
		self.list = [None, None, None]

	def __init__(self, x, y, z):
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
		self.listc = [x, y, z]

	def __repr__(self):
		s = "Point(%s, %s, %s)" % (self.x, self.y, self.z)
		return s	

	def __str__(self):
		s = "Point(%s, %s, %s)" % (self.x, self.y, self.z)
		return s

	def __mul__(self, num):
		return Point(self.x*num, self.y*num, self.z*num)

	def __add__(self, other):
		return Point(self.x+other.x, self.y+other.y, self.z+other.z)

	def __sub__(self,other):
		return Point(self.x-other.x, self.y-other.y, self.z-other.z)

	def __truediv__(self, num):
		return Point(self.x/num, self.y/num, self.z/num)

	def __div__(self, num):
		return Point(self.x/num, self.y/num, self.z/num)

	def __getitem__(self, key):                  
		if self.listc is None:
			raise TypeError('not indexable')
		return self.listc[key]

	def distance(self, other):
		return np.sqrt((self.x-other.x)**2+(self.y-other.y)**2+(self.z-other.z)**2)



