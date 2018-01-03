from math import *
import math
from cmath import exp as Cexp
import scipy.linalg as linalg
import numpy
from scipy import integrate


#################################################
#########  fundamental constants  ###############
hbar = 0.65821e-15#Dirac's constant in eV . s
h = 4.13566743e-15#Planck's constant in eV.s
e = 1.60217662e-19 #Quantum of charge in C
kB = 0.00008617 #Boltzmann's constant in eV/K
################################################
################################################

###################################################
####### geometry parameters phosphorene ###########
a1 = 2.22# Angstroms
a2 = 2.24# Angstroms
alpha1 = math.radians(96.5)
alpha2 = math.radians(101.9)
beta = math.radians(72)
d1 = a1 * sin(alpha1/2)
d2 = a1 * cos(alpha1/2) + a2 * cos(beta) 
d3 = a1 * cos(alpha1/2)
d4 = a2 * cos(beta)
l1 = 2 * d1
l2 = 2 * d2
###################################################
################################################

##################################
####  hoppings in eV  ############
t1 = -1.486
t2 = 3.729
t3 = -0.252
t4 = -0.071
t5 = -0.019
t6 = 0.186
t7 = -0.063
t8 = 0.101
t9 = -0.042
t10 = 0.073
##################################
##################################


######################################## structure factors #########################################################################################	
def tAA(kx, ky):
	value = 2 * t3 * cos(2 * d1 * kx) + 2 * t7 * cos(2 * d2 * ky) + 4 * t10 * cos(2 * d1 * kx) * cos(2 * d2 * ky)
	return value

def tAB(kx, ky):
	value = 2 * t1 * cos(d1 * kx) * Cexp(-1j * d3 * ky) + 2 * t4 * cos(d1 * kx) * Cexp(1j * (d2 + d4) * ky) + 2 * t8 * cos(3 * d1 * kx) * Cexp(-1j * d3 * ky)
	return value

def tAC(kx, ky):
	value = t2 * Cexp(1j * d4 * ky) + t6 * Cexp(-1j * (d2 + d3) * ky) + 2 * t9 * cos(2 * d1 * kx) * Cexp(-1j * (d2 + d3) * ky)
	return value

def tAD(kx, ky):
	value = 4 * t5 * cos(d1 * kx) * cos(d2 * ky)
	return value
########################################################################################################################################################

################# Tight-binding Hamiltonian, eigenvectors and eigenvalues ##############################################################################
def H(kx, ky):
	hamiltonian = numpy.matrix([[tAA(kx, ky) + tAD(kx, ky), tAB(kx, ky).conjugate() + tAC(kx, ky).conjugate()],[tAB(kx, ky) + tAC(kx, ky) , tAA(kx, ky) + tAD(kx, ky)]])
	return hamiltonian	
# eigenvalues
def E0(kx, ky):
	return linalg.eigh(H(kx, ky))[0][0]

def E1(kx, ky):
	return linalg.eigh(H(kx, ky))[0][1]

# normalized eigenvector correponding to the eigenvalue n : n = 0 or 1, for E0 or E1, respectively
def Psi(kx, ky,n):
	eigenvector = linalg.eigh(H(kx, ky))
	return eigenvector[1][:,n]
#########################################################################################################################################################


#Here, I define the amplitude F =  \<psi k, n | psi k + q, m >| **2 which goes in the integrand of the pair-bubble integral
#qx and qy are the components of the momentum transfer, whereas n and m defines the energy bands
#The only restriction is that it is not normalized in this first version of the code.
def F(kx, ky, qx, qy, n, m):
	return numpy.dot(Psi(kx, ky, n).conjugate(), Psi(kx + qx, ky + qy, m))


########## Fermi-Dirac distribution #######################
# T is the temperature given in K
# EF is the chemical potential
def f0(kx, ky, T, Ef):
	b = 1/(kB * T)
	fermi = 1/(exp(b * (E0(kx, ky) - E0(0, 0) - Ef)) + 1)
	return fermi

def f1(kx, ky, T, Ef):
	b = 1/(kB * T)
	fermi = 1/(exp(b * (E1(kx, ky) - E1(0, 0) - Ef)) + 1)
	return fermi
##########################################################

############ Argument of the Pair-Bubble ###############################
def Argument(kx, ky, qx, qy, w, T, Ef, n, m):
	pairbubble = - 2/(2 * pi)**2 * (f0(kx, ky, T, Ef) - f0(kx + qx, ky + qy, T, Ef))/(E0(kx, ky) - E0(kx + qx, ky + qy) + hbar * w) * F(kx, ky, qx, qy, n, m)
	return pairbubble

#def Pair_Bubble(qx, qy, w):
#	return integrate.nquad(Argument, [[-inf, +inf],[-inf, +inf]])[0]

#def Dieletric_function(qx, qy, w):
#	return 1 + e**2/(2 * sqrt(qx**2 + qy**2)) * Pair_Bubble(qx, qy, w)

Qx = [ 0.2 * i/100 for i in range(100)]
Qy = [ 0.2 * i/100 for i in range(100)]
Kx = [-2 + 1e-100 + 4 * i/500 for i in range(500)]
Ky = [-2 + 1e-100 + 4 * i/500 for i in range(500)]

file1 = open('dieletric_functiony_TB.dat','w')
#for kx in xmomentum:
#	for ky in ymomentum:
#		E = E1(kx, ky) - E1(0, 0)
#		f = f1(kx, ky, 200, 0)
#		out = '%e    %e\n' % (E, f)
#		file1.write(out)
#		print(ky)
##########################################################################
#for qx in Qx:
#	for qy in Qy:
#		print(1 - qx/2)
#		w = 1e-100
#		T = 300
#		q =  sqrt(qx**2 + qy**2)		
#		def argument(kx, ky):
#			return Argument(kx, ky, qx, qy, w, T, Ef = 0, n = 1, m = 1)
#		def arg(kx):
#			return sum([argument(kx, ky) * 2/500 for ky in Ky])
#		pair_bubble = sum([arg(kx) * 2/500 for kx in Kx])
#		dieletric_function = 1 + e**2/(2 * epsilon * q) * pair_bubble		
#		out = '%e    %e\n' % (q, dieletric_function)
#		file1.write(out)
#file1.close()
############################################################################
for qy in Qx:
	qx = 0
	print(1 - qy/2)
	w = 1e-100
	T = 300
	q =  sqrt(qx**2 + qy**2)		
	def argument(kx, ky):
		return Argument(kx, ky, qx, qy, w, T, Ef = 0, n = 1, m = 1)
	def arg(kx):
		return sum([argument(kx, ky) * 2/500 for ky in Ky])
	pair_bubble = sum([arg(kx) * 2/500 for kx in Kx])
	dieletric_function = 1/(2 * (q + 1e-100)) * pair_bubble		
	out = '%e    %e\n' % (q, dieletric_function)
	file1.write(out)
file1.close()
		
