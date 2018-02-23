import mol
import numpy as np
import integrals
import scipy.linalg as slin

# SCF Class which houses the SCF procedure and all of the matrix building functions
# Once the scf object runs the kernel, the user has access to any of the matrices and energies
# This is currently the highest level in the folder, which the driver file will directly call
class SCF:

	def __init__(self, mol):
		self.mol = mol
		self.S = None
		self.C = None
		self.F = None
		self.e = None
		self.P = None
		self.T = None
		self.Hcore = None
		self.eris = None
		self.Vne = None
		self.Vee = None
		self.e_e = None
		self.z_e = None
		self.e_tot = None

# The kernel function for the SCF class, this is the meat of the SCF procedure
# Here all of the matrices are built and plugged into the SCF routine
# The end result will be adding all pertinent matrices to the SCF object and
# calculating the energy of the system
	def kernel(self):
		nmo = self.mol.nmo()
		self.S = self.create_S()
		self.T = self.create_T()
		self.Vne = self.create_enuc()
		self.eris = self.create_eris()
		self.Hcore = self.T + self.Vne
		self.P = np.zeros((nmo, nmo))
		P_old = self.P
	
		diff = False
		while not diff:
			#self.Vee = self.create_Vee(self.mol, eris)#, P)
			self.Vee = self.create_Vee()
			self.F = self.Hcore + self.Vee
	
			self.e, self.C = slin.eig(self.F, self.S)
			idx = self.e.argsort()
			self.e = self.e[idx]
			self.C = self.C[:, idx]
			Nf = np.zeros(len(self.C))
			for m in range(len(self.C)):
				for i in range(len(self.C)):
					for j in range(len(self.C)):
							Nf[m] += self.C[i, m]*self.C[j, m]*self.S[i, j]
			for m in range(len(self.C)):
				self.C[:, m] = self.C[:, m]/np.sqrt(Nf[m])
	
			self.P = self.create_P()
			diff = np.allclose(self.P, P_old, atol=1e-6)
			P_old = self.P
		
		self.e_e = self.e_energy()
		self.z_e = self.Z_energy()
		self.e_tot = self.tot_energy() 
	
	# Function to build the density matrix (P) from the C matrix
	# Input: self (uses self.C)
	# Output: out_P -> density matrix from current C
	def create_P(self):
		nmo = self.mol.nmo()
		out_P = np.zeros((nmo, nmo))
		mos = []
		elec = self.mol.elec
		for atom in self.mol.atoms:
			for orb in atom.orbitals:
				mos.append(orb)
		for muidx in range(nmo):
			for nuidx in range(nmo):
				sum_add = 0
				for a in range(int(elec/2)):
					sum_add += 2*self.C[muidx, a]*self.C[nuidx, a]
				out_P[muidx, nuidx] = sum_add
		return out_P
	# Function to build the overlap matrix (S) from scratch
	# Input: self
	# Output out_S -> overlap matrix for molecule in system
	def create_S(self):
		nmo = self.mol.nmo()
		out_S = np.zeros((nmo, nmo))
		mos = []
		for atom in self.mol.atoms:
			for orb in atom.orbitals:
				mos.append(orb)
		for muidx in range(len(mos)):
			mu = mos[muidx]
			for nuidx in range(muidx, len(mos)):
				nu = mos[nuidx]
				out_S[muidx, nuidx] = out_S[nuidx, muidx] = integrals.overlap(mu, nu)
		return out_S

	# Function to build the kinetic energy matrix (T) from scratch
	# Input: self
	# Output: out_T -> kinetic energy matrix built from scratch
	def create_T(self):
		nmo = self.mol.nmo()
		out_T = np.zeros((nmo, nmo))
		mos = []
		for atom in self.mol.atoms:
			for orb in atom.orbitals:
				mos.append(orb)
		for muidx in range(len(mos)):
			mu = mos[muidx]
			for nuidx in range(muidx, len(mos)):
				nu = mos[nuidx]
				out_T[muidx, nuidx] = out_T[nuidx, muidx] = integrals.kinetic(mu, nu)
		return out_T

	# Function to build the electron nuclear attraction matrix (Vne) from scratch
	# Input: self
	# Output: out_enuc -> electron nuclear attraction matrix
	def create_enuc(self):
		nmo = self.mol.nmo()
		out_enuc = np.zeros((nmo, nmo))
		mos = []
		for atom in self.mol.atoms:
			for orb in atom.orbitals:
				mos.append(orb)
		for atom in self.mol.atoms:
			for muidx in range(len(mos)):
				mu = mos[muidx]
				for nuidx in range(muidx, len(mos)):
					nu = mos[nuidx]
					to_sub = atom.Z*integrals.elec_nuc(mu, nu, atom)
					out_enuc[muidx, nuidx] -= to_sub
					if nuidx != muidx:
						out_enuc[nuidx, muidx] -= to_sub
		return out_enuc
		
	# Function to build the four-indexed ERIS (NOTE: This isn't the same as the electron-electron repulsion matrix, Vee)
	# Input: self
	# Output: out_eris -> four-indexed electron repulsion integrals (eris)
	def create_eris(self):
		nmo = self.mol.nmo()
		out_eris = np.zeros((nmo, nmo, nmo, nmo))
		mos = []
		for atom in self.mol.atoms:
			for orb in atom.orbitals:
				mos.append(orb)
		for muidx in range(len(mos)):
			mu = mos[muidx]
			for nuidx in range(muidx, len(mos)):
				nu = mos[nuidx]
				for kapidx in range(muidx, len(mos)):
					kap = mos[kapidx]
					for lamidx in range(kapidx, len(mos)):
						lam = mos[lamidx]
						out_eris[muidx, nuidx, kapidx, lamidx] = out_eris[nuidx, muidx, kapidx, lamidx] = out_eris[muidx, nuidx, lamidx, kapidx] = out_eris[nuidx, muidx, lamidx, kapidx] = out_eris[kapidx, lamidx, muidx, nuidx] = out_eris[kapidx, lamidx, nuidx, muidx] = out_eris[lamidx, kapidx, muidx, nuidx] = out_eris[lamidx, kapidx, nuidx, muidx] = integrals.eris(mu, nu, kap, lam)
		return out_eris
	
	# Function to build the electron-electron repulsion matrix (Vee) from the eris
	# Input self (uses self.eris)
	# Output: out_vee -> electron-electron repulsion matrix (Vee)
	def create_Vee(self):
		nmo = self.mol.nmo()
		out_vee = np.zeros((nmo, nmo))
		mos = []
		for atom in self.mol.atoms:
			for orb in atom.orbitals:
				mos.append(orb)
		for muidx in range(len(mos)):
			for nuidx in range(len(mos)):
				sum_add = 0
				for kapidx in range(len(mos)):
					for lamidx in range(len(mos)):
						sum_add += self.P[kapidx, lamidx]*(self.eris[muidx, nuidx, lamidx, kapidx] - 0.5*self.eris[muidx, kapidx, lamidx, nuidx])
				out_vee[muidx, nuidx] = sum_add
		return out_vee					
	
	# Function to determine the electronic energy (electronic part of Hamiltonian)
	# Input: self (uses self.P, self.Hcore, self.F)
	# Output: sum_add -> HF electronic energy 
	def e_energy(self):
		sum_add = 0
		for mu in range(len(self.P[0])):
			for nu in range(len(self.P[0])):
				sum_add += self.P[mu, nu]*(self.Hcore[mu, nu] + self.F[mu, nu])
		return 0.5*sum_add

	# Function to determine the nuclear-nuclear repulsion energy
	# Input: self
	# Output: sum_add -> nuclear-nuclear repulsion energy
	def Z_energy(self):
		sum_add = 0
		for A in self.mol.atoms:
			for B in self.mol.atoms:
				if A != B:
					sum_add += A.Z*B.Z/A.distance(B)
		return 0.5*sum_add

	# Function which adds the electronic and nuclear parts of the energy for the total energy
	# Input: self (uses self.e_e, self.z_e)
	# Output: total energy of system
	def tot_energy(self):
		return self.e_e + self.z_e

