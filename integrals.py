import numpy as np
import scipy
import copy
import sys
import scipy.misc as smisc
# Overlap matrix element parent function which calls required primitive integral functions to determine S_(mu,nu)
# This function should be called as one builds the overlap (S) matrix
# Input: mu,nu -> atomic orbital objects
# Output: sum_add -> overlap matrix element, S_(mu,nu)
def overlap(mu, nu):
	sum_add = 0
	for k in mu.pgo:
		for l in nu.pgo:
			a = mu.l#[mu.l, mu.l, mu.l]
			b = nu.l#[nu.l, nu.l, nu.l]
			alpha = k.exp
			beta = l.exp
			A = mu.xyz
			B = nu.xyz
			P = P_gen(alpha, beta, A, B)
			KAB = kab_gen(alpha, beta, A, B)
			p = p_gen(alpha, beta)
			sum_add += k.coeff*l.coeff*norm(a, alpha)*norm(b, beta)*prim_3D_overlap(a, b, P, A, B, KAB, p)
	return sum_add

# Kinetic energy matrix element parent function which calls required primitive integral functions to determine K_(mu,nu)
# This function should be called as one builds the kinteic energy (T) matrix
# Input: mu,nu -> atomic orbital objects
# Output: sum_add -> kinetic energy matrix element, T_(mu,nu)
def kinetic(mu, nu):
	sum_add = 0
	for k in mu.pgo:
		for l in nu.pgo:
			a = mu.l#[mu.l, mu.l, mu.l]
			b = nu.l#[nu.l, nu.l, nu.l]
			alpha = k.exp
			beta = l.exp
			A = mu.xyz
			B = nu.xyz
			P = P_gen(alpha, beta, A, B)
			KAB = kab_gen(alpha, beta, A, B)
			p = p_gen(alpha, beta)
			sum_add += k.coeff*l.coeff*norm(a, alpha)*norm(b, beta)*prim_3D_kinetic(a, b, P, A, B, KAB, p, beta)
	return sum_add

# Electron-nuclear attraction matrix element parent function which calls required primitive integral functions to determine K_(mu,nu)
# This function should be called as one builds the electron-nuclear attraction (Enuc) matrix
# Input: mu,nu -> atomic orbital objects; nuc -> nucleus coordinates
# Output: sum_add -> electron-nuclear attraction matrix element, Enuc_(mu,nu)
def elec_nuc(mu, nu, nuc):
	sum_add = 0
	for k in mu.pgo:
		for l in nu.pgo:
			a = mu.l
			b = nu.l
			alpha = k.exp
			beta = l.exp
			A = mu.xyz
			B = nu.xyz
			p = p_gen(alpha, beta)
			P = P_gen(alpha, beta, A, B)
			KAB = kab_gen(alpha, beta, A, B)
			C = nuc.xyz
			m = 0
			sum_add += k.coeff*l.coeff*norm(a, alpha)*norm(b, beta)*prim_enuc(a, b, A, B, C, P, KAB, p, m)
	return sum_add

# Function which calculates the four-indexed electron repulsion integral (mu,nu|kappa,lambda)
# This function should be called as one builds the eris and electron-electron repulsion (Vee) matrix
# Input: mu,nu,kap,lam -> atomic orbital objects
# Output: sum_add -> four-indexed electron repulsion integral, (mu,nu|kap,lam)
def eris(mu, nu, kap, lam):
	sum_add = 0
	for k in mu.pgo:
		for l in nu.pgo:
			for n in kap.pgo:
				for o in lam.pgo:
					a = mu.l#[mu.l, mu.l, mu.l]
					b = nu.l#[nu.l, nu.l, nu.l]
					c = kap.l#[kap.l, kap.l, kap.l]
					d = lam.l#[lam.l, lam.l, lam.l]
					A = mu.xyz
					B = nu.xyz
					C = kap.xyz
					D = lam.xyz	
					alpha = k.exp
					beta = l.exp
					gamma = n.exp
					delta = o.exp
					p = p_gen(alpha, beta)
					q = p_gen(gamma, delta)
					P = P_gen(alpha, beta, A, B)
					Q = P_gen(gamma, delta, C, D)
					W = P_gen(p, q, P, Q)
					KAB = kab_gen(alpha, beta, A, B)
					KCD = kab_gen(gamma, delta, C, D)
					m = 0
					sum_add += k.coeff*l.coeff*n.coeff*o.coeff*norm(a, alpha)*norm(b, beta)*norm(c, gamma)*norm(d, delta)*prim_eris(a, b, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m)
	return sum_add

# 3D primitive integrals used for the overlap and kinetic energy matrix elements
# Input: a,b -> PGO angular momentum; A,B -> PGO coordinates; P -> weighted midpoint coordinates; KAB -> Gaussian pre-factor; p -> weighted exponent
# Output: 3D primitive integral for given pair of PGOs
def prim_3D_overlap(a, b, P, A, B, KAB, p):
	return (np.pi/p)**(3/2)*KAB*prim_1D_overlap(a[0], b[0], P.x, A.x, B.x, p)*prim_1D_overlap(a[1], b[1], P.y, A.y, B.y, p)*prim_1D_overlap(a[2], b[2], P.z, A.z, B.z, p)

# 1D primitive integrals used for the overlap and kinetic energy matrix elements
# Input: aw,bw -> PGO angular momentum (w coord); Aw,Bw -> PGO coordinates (w coord); Pw -> weighted midpoint coordinates (w coord); p -> weighted exponent
# Output: 1D (w dimension) primitive integral for given pair of PGOs
def prim_1D_overlap(aw, bw, Pw, Aw, Bw, p):
	sum_add = 0
	i = 0
	while i <= (aw+bw)/2:
		sum_add += fk(2*i, aw, bw, Pw, Aw, Bw)*df(2*i-1)/(2*p)**i
		i += 1
	return sum_add
# Permutation function for the 1D integrals
# Input: kidx -> index for summation bounds; aw,bw -> PGO angular momentum (w coord); Aw,Bw -> PGO coordinates Pw -> weighted midpoint coordinates (w coord)
# Output: sum_add -> permutation function output
def fk(kidx, aw, bw, Pw, Aw, Bw):
	sum_add = 0
	j = max(0, kidx-aw)
	while j <= min(kidx, bw):
		sum_add += smisc.comb(aw, kidx-j)*smisc.comb(bw, j)*(Pw-Aw)**(aw-kidx+j)*(Pw-Bw)**(bw-j)
		j += 1
	return sum_add

# Normalization function for given PGO
# Input: a -> PGO angular momentum, alpha -> PGO exponent
# Output: Normalization coefficient for given PGO
def norm(a, alpha):
	return (2/np.pi)**(3/4)*2**(sum(a))*alpha**((2*sum(a)+3)/4)/np.sqrt(df(2*a[0]-1)*df(2*a[1]-1)*df(2*a[2]-1))

# Double factorial (x)!! function
# Input: x -> integer
# Output: out -> double factorial of x, x!!
def df(x):
	out = 1
	while x > 0:
		out *= x
		x -= 2
	return out

# Function to calculate the Gaussian product prefactor, KAB
# Input: alpha,beta -> PGO exponents; A,B -> PGO coordinates
# Output: Gaussian product prefactor, KAB
def kab_gen(alpha, beta, A, B):
	return np.exp((-alpha*beta/(alpha+beta))*A.distance(B)**2)

# Exponent sum for PGO product
# Input: alpha,beta -> PGO exponenets
# Output: sum of PGO exponents
def p_gen(alpha, beta):
	return alpha+beta

# Weighted midpoint of coordinates for PGO pair
# Input: alpha,beta -> PGO exponents; A,B -> PGO coordinates
# Output: weighted midpoint for PGO pair (P)
def P_gen(alpha, beta, A, B):
	return (A*alpha + B*beta)/(alpha+beta)

# 3D primitive integrals for kinetic energy matrix
# Input: a,b -> PGO angular momentum; A,B -> PGO coordinates; P -> weighted midpoint coordinates; KAB -> Gaussian pre-factor; p -> weighted exponent; beta -> PGO exponent
# Output: 3D primitive integral for given pair of PGOs for kinetic energy matrix
def prim_3D_kinetic(a, b, P, A, B, KAB, p, beta):
	return prim_1D_kinetic(a, b, P, A, B, KAB, p, beta, 0) + prim_1D_kinetic(a, b, P, A, B, KAB, p, beta, 1) + prim_1D_kinetic(a, b, P, A, B, KAB, p, beta, 2)


# 1D primitive integrals used for the overlap and kinetic energy matrix elements, which calls 3D primitive integrals for overlap matrix
# Input: a,b -> PGO angular momentum; A,B -> PGO coordinates; P -> weighted midpoint coordinates; KAB -> Gaussian pre-factor; p -> weighted exponent; beta -> PGO exponent
# Input: w -> dimension of interest (0,1,2)
# Output: 1D primitive integral for given pair of PGOs for kinetic energy matrix
def prim_1D_kinetic(a, b, P, A, B, KAB, p, beta, w):
	bPlus = copy.copy(b)
	bPlus[w] += 2
	bMinus = copy.copy(b)
	bMinus[w] -= 2
	return beta*(2*b[w]+1)*prim_3D_overlap(a, b, P, A, B, KAB, p) - 2*beta**2*prim_3D_overlap(a, bPlus, P, A, B, KAB, p) - 0.5*b[w]*(b[w]-1)*prim_3D_overlap(a, bMinus, P, A, B, KAB, p)

# Recursion parent function for the enuc integrals
# Input: a,b -> PGO angular momentum; A,B -> PGO coordinates; P -> weighted midpoint coordinates; KAB -> Gaussian pre-factor; p -> weighted exponent; m -> order of auxiliary integral
# Output: enux integral element for given PGOs, which added together will give Enuc_(mu,nu)
def prim_enuc(a, b, A, B, C, P, KAB, p, m):
	if a == [0,0,0] and b == [0,0,0]:
		T = p*(P.distance(C))**2
		return 2*np.pi/p*KAB*boys(m,T)
	elif any(n < 0 for n in a) or any(n < 0 for n in b):
		return 0
	elif a[0] != 0:
		return enuc_single('a', 0, a, b, A, B, C, P, KAB, p, m)
	elif a[1] != 0:
		return enuc_single('a', 1, a, b, A, B, C, P, KAB, p, m)
	elif a[2] != 0:
		return enuc_single('a', 2, a, b, A, B, C, P, KAB, p, m)
	elif b[0] != 0:
		return enuc_single('b', 0, a, b, A, B, C, P, KAB, p, m)
	elif b[1] != 0:
		return enuc_single('b', 1, a, b, A, B, C, P, KAB, p, m)
	elif b[2] != 0:
		return enuc_single('b', 2, a, b, A, B, C, P, KAB, p, m)
	
# Enuc integral single function used for reducing both a and b angular momenta and all three dimensions
# Input: a,b -> PGO angular momentum; A,B -> PGO coordinates; P -> weighted midpoint coordinates; KAB -> Gaussian pre-factor; p -> weighted exponent; m -> order of auxiliary integral
# Input: g -> angular momentum of interest (a,b), w -> dimension of interest (0,1,2)
# Output: recursion integral for enuc recursion relation for given dimension and angular momentum 
def enuc_single(g, w, a, b, A, B, C, P, KAB, p, m):
	aw = a[w]
	bw = b[w]
	aMinus = copy.copy(a)
	aMinus[w] -= 1
	bMinus = copy.copy(b)
	bMinus[w] -= 1
	if g == 'a':
		aDminus = copy.copy(a)
		aDminus[w] -= 2
		return (P-A)[w]*prim_enuc(aMinus, b, A, B, C, P, KAB, p, m)+(C-P)[w]*prim_enuc(aMinus, b, A, B, C, P, KAB, p, m+1)+(aw-1)/(2*p)*(prim_enuc(aDminus, b, A, B, C, P, KAB, p, m)-prim_enuc(aDminus, b, A, B, C, P, KAB, p, m+1))+bw/(2*p)*(prim_enuc(aMinus, bMinus, A, B, C, P, KAB, p, m) - prim_enuc(aMinus, bMinus, A, B, C, P, KAB, p, m+1))
	else:
		bDminus = copy.copy(b)
		bDminus[w] -= 2
		return (P-B)[w]*prim_enuc(a, bMinus, A, B, C, P, KAB, p, m)+(C-P)[w]*prim_enuc(a, bMinus, A, B, C, P, KAB, p, m+1)+(aw-1)/(2*p)*(prim_enuc(a, bDminus, A, B, C, P, KAB, p, m)-prim_enuc(a, bDminus, A, B, C, P, KAB, p, m+1))+bw/(2*p)*(prim_enuc(aMinus, bMinus, A, B, C, P, KAB, p, m) - prim_enuc(aMinus, bMinus, A, B, C, P, KAB, p, m+1))
	
# Boys function used in enuc and eris integrals
# Input: m -> order of auxiliary integral, T -> weighted weighted midpoint
# Output: Boys function with given m and T
def boys(m, T):
#	if T == 0:
	if T < 1e-16:
		return 1/(2*m+1)
	else:
		return scipy.special.gamma(m+0.5)*scipy.special.gammainc(m+0.5, T)/(2*T**(m+0.5))

# Electron repulsion integral recursion relation parent function which calls horizontal relations
# Input: a,b,c,d -> PGO angular momentum; A,B,C,D -> PGO coordinates; P,Q,W -> weighted midpoint coordinates; KAB,KCD -> Gaussian pre-factor; p,q -> weighted exponent; m -> order of auxiliary integral
# Output: eri for given combination of four PGOs
def prim_eris(a, b, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m):
	if a == [0,0,0] and b == [0,0,0] and c == [0,0,0] and d == [0,0,0]:
		T = p*q/(p+q)*P.distance(Q)**2
		return 2*np.pi**(5/2)/(p*q*np.sqrt(p+q))*KAB*KCD*boys(m, T)
	elif any(n < 0 for n in a) or any(n < 0 for n in b) or any(n < 0 for n in c) or any(n < 0 for n in d):
		return 0
	else:
		return prim_eris_hrr(a, b, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m)

# Vertical eris recursion parent function which displaces a,c angular momentum to auxiliary integral order, m
# Input: a,c -> PGO angular momentum; A,C -> PGO coordinates; P,Q,W -> weighted midpoint coordinates; KAB,KCD -> Gaussian pre-factor; p,q -> weighted exponent; m -> order of auxiliary integral
# Output: vertical recursion relation for eris
def prim_eris_vrr(a, c, A, C, p, q, P, Q, W, KAB, KCD, m):
	if a == [0,0,0] and c == [0,0,0]:
		T = p*q/(p+q)*P.distance(Q)**2
		return 2*np.pi**(5/2)/(p*q*np.sqrt(p+q))*KAB*KCD*boys(m, T)
	elif any(n < 0 for n in a) or any(n < 0 for n in c):
		return 0
	elif a[0] != 0:
		return eris_vrr_single('a', 0, a, c, A, C, p, q, P, Q, W, KAB, KCD, m)	
	elif a[1] != 0:
		return eris_vrr_single('a', 1, a, c, A, C, p, q, P, Q, W, KAB, KCD, m)	
	elif a[2] != 0:
		return eris_vrr_single('a', 2, a, c, A, C, p, q, P, Q, W, KAB, KCD, m)	
	elif c[0] != 0:
		return eris_vrr_single('c', 0, a, c, A, C, p, q, P, Q, W, KAB, KCD, m)
	elif c[1] != 0:
		return eris_vrr_single('c', 1, a, c, A, C, p, q, P, Q, W, KAB, KCD, m)
	elif c[2] != 0:
		return eris_vrr_single('c', 2, a, c, A, C, p, q, P, Q, W, KAB, KCD, m)

# Vertical recursion relation for given dimension and chosen a,c angular momentum		
# Input: a,c -> PGO angular momentum; A,C -> PGO coordinates; P,Q,W -> weighted midpoint coordinates; KAB,KCD -> Gaussian pre-factor; p,q -> weighted exponent; m -> order of auxiliary integral
# Input: g -> angular momentum of interest (a,b), w -> dimension of interest (0,1,2)
# Output: vertical recursion relation for given a,c angular momentum, dimension, w, and auxiliary order, m
def eris_vrr_single(g, w, a, c, A, C, p, q, P, Q, W, KAB, KCD, m):	
	aMinus = copy.copy(a)
	cMinus = copy.copy(c)
	aMinus[w] -= 1
	cMinus[w] -= 1
	if g == 'a':
		aDminus = copy.copy(a)
		aDminus[w] -= 2	
		return (P-A)[w]*prim_eris_vrr(aMinus, c, A, C, p, q, P, Q, W, KAB, KCD, m) + (W-P)[w]*prim_eris_vrr(aMinus, c, A, C, p, q, P, Q, W, KAB, KCD, m+1) + (a[w]-1)/(2*p)*(prim_eris_vrr(aDminus, c, A, C, p, q, P, Q, W, KAB, KCD, m) - q/(p+q)*prim_eris_vrr(aDminus, c, A, C, p, q, P, Q, W, KAB, KCD, m+1)) + c[w]/(2*(p+q))*prim_eris_vrr(aMinus, cMinus, A, C, p, q, P, Q, W, KAB, KCD, m+1)
	else:
		cDminus = copy.copy(c)
		cDminus[w] -= 2
		return (Q-C)[w]*prim_eris_vrr(a, cMinus, A, C, p, q, P, Q, W, KAB, KCD, m) + (W-Q)[w]*prim_eris_vrr(a, cMinus, A, C, p, q, P, Q, W, KAB, KCD, m+1) + (c[w]-1)/(2*q)*(prim_eris_vrr(a, cDminus, A, C, p, q, P, Q, W, KAB, KCD, m) - p/(p+q)*prim_eris_vrr(a, cDminus, A, C, p, q, P, Q, W, KAB, KCD, m+1)) + a[w]/(2*(p+q))*prim_eris_vrr(aMinus, cMinus, A, C, p, q, P, Q, W, KAB, KCD, m+1)

# Horizontal recursion relation for given dimension and chosen b,d angular momentum		
# Input: a,b,c,d -> PGO angular momentum; A,B,C,D -> PGO coordinates; P,Q,W -> weighted midpoint coordinates; KAB,KCD -> Gaussian pre-factor; p,q -> weighted exponent; m -> order of auxiliary integral
# Input: g -> angular momentum of interest (a,b), w -> dimension of interest (0,1,2)
# Output: Horizontal recursion relation for eris for given b,d, dimension, w, and auxiliary integral order, m
def eris_hrr_single(g, w, a, b, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m):
	aPlus = copy.copy(a)
	bMinus = copy.copy(b)
	cPlus = copy.copy(c)
	dMinus = copy.copy(d)
	aPlus[w] += 1
	bMinus[w] -= 1
	cPlus[w] += 1
	dMinus[w] -= 1
	if g == 'b':
		return prim_eris_hrr(aPlus, bMinus, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m) + (A-B)[w]*prim_eris_hrr(a, bMinus, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m)
	else:
		return prim_eris_hrr(a, b, cPlus, dMinus, A, B, C, D, p, q, P, Q, W, KAB, KCD, m) + (C-D)[w]*prim_eris_hrr(a, b, c, dMinus, A, B, C, D, p, q, P, Q, W, KAB, KCD, m)
		
# Horizontal eris recursion parent function which displaces b,d angular momentum to a,c angular momentum
# Input: a,b,c,d -> PGO angular momentum; A,B,C,D -> PGO coordinates; P,Q,W -> weighted midpoint coordinates; KAB,KCD -> Gaussian pre-factor; p,q -> weighted exponent; m -> order of auxiliary integral
# Output: vertical recursion relation for eris
def prim_eris_hrr(a, b, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m):
	if a == [0,0,0] and b == [0,0,0] and c == [0,0,0] and d == [0,0,0]:
		T = p*q/(p+q)*P.distance(Q)**2
		return 2*np.pi**(5/2)/(p*q*np.sqrt(p+q))*KAB*KCD*boys(m, T)
	elif b == [0,0,0] and d == [0,0,0]:
		return prim_eris_vrr(a, c, A, C, p, q, P, Q, W, KAB, KCD, m)
	elif any(n < 0 for n in b) or any(n < 0 for n in d):
		return 0
	elif b[0] != 0:
		return eris_hrr_single('b', 0, a, b, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m)
	elif b[1] != 0:
		return eris_hrr_single('b', 1, a, b, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m)
	elif b[2] != 0:
		return eris_hrr_single('b', 2, a, b, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m)
	elif d[0] != 0:
		return eris_hrr_single('d', 0, a, b, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m)
	elif d[1] != 0:
		return eris_hrr_single('d', 1, a, b, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m)
	elif d[2] != 0:
		return eris_hrr_single('d', 2, a, b, c, d, A, B, C, D, p, q, P, Q, W, KAB, KCD, m)


