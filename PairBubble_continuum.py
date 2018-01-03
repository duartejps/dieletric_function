from math import *
import math
from cmath import exp as Cexp
import scipy.linalg as linalg
import numpy
from scipy import integrate

#########################################
#In this code, the fundamental units are:
#Energy : eV
#distance : Angstrom
#Time : second
#Vacuum permittivity : 1 / 4 * pi


#########  fundamental constants  ###############
hbar = 0.65821e-15#Dirac's constant in eV . s
h = 4.13566743e-15#Planck's constant in eV.s
e = 1.60217662e-19 #Quantum of charge in C
kB = 0.00008617 #Boltzmann's constant in eV/K
c = 299792458#velocity of light in m/s
alpha = 1/137.036#fine-structure constant
################################################

#-----------
dAA = -0.338
dAB = -2.912
dAC = 3.831
dAD = -0.076

nAA = 1.161
nAB = 2.05
nAC = 0.460
nAD = 0.104

yAA = -1.563
yAB = 3.607
yAC = -1.572
yAD = 0.179

xAB = 3.688
xAC = 2.208
#----------
#--------------
u0 = dAA + dAD
nx = nAA + nAD
yx = nAB + nAC
d = dAB + dAC
ny = yAA + yAD
yy = yAB + yAC
x = xAB + xAC
#-------------

# Energy bands calculated from the continuum approximation. n = +1 for electrons and -1 for holes.
def E(kx, ky, n):
	return u0 + nx * kx**2 + ny * ky**2 + n * sqrt((d + yx * kx**2 + yy * ky**2)**2 + (x * ky)**2)
#some functions used to calculate the eigenspinor associated with the energy eigenvalues given above.
def fv(kx, ky):
	return u0 - d + (nx - yx) * kx**2 + (ny - yy) * ky**2
def fc(kx, ky):
	return u0 + d + (nx + yx) * kx**2 + (ny + yy) * ky**2
def theta(kx, ky):
	return atan(2 * x * ky / (fc(kx, ky) - fv(kx, ky)))
def F(kx, ky, qx, qy, n, m):
	return 0.5 * (1 + n * m * cos(theta(kx + qx, ky + qy) - theta(kx, ky)))
########## Fermi-Dirac distribution #######################
# T is the temperature given in K
# EF is the chemical potential
def f(kx, ky, T, n, Ef):
	b = 1/(kB * T)
	try:
		fermi = 1/(exp(b * (E(kx, ky, n) - E(0, 0, n) - Ef)) + 1)
	except OverflowError:
		fermi = float('inf')	
	return fermi
##########################################################

############ Argument of the Pair-Bubble ###############################
def Argument_real(kx, ky, qx, qy, w, T, Ef, n, m):
	pairbubble =  (2/(4 * pi ** 2)) * (f(kx + qx, ky + qy, T, n, Ef) - f(kx, ky, T, n, Ef))/(E(kx, ky, n) - E(kx + qx, ky + qy, n) + hbar * w) * F(kx, ky, qx, qy, n, m)
	return pairbubble

def delta_dirac(kx, ky, qx, qy, n, w):
	epsilon = 1e-100
	delta = 1/(sqrt(4 * pi* epsilon)) * exp(-((E(kx, ky, n) - E(kx + qx, ky + qy, n) + hbar * w)**2)/(4 * epsilon))

def Argument_imag(kx, ky, qx, qy, w, T, Ef, n, m):
	pairbubble =  (1/(2 * pi)) * (f(kx + qx, ky + qy, T, n, Ef) - f(kx, ky, T, n, Ef)) * delta_dirac(kx, ky, qx, qy, n, w) * F(kx, ky, qx, qy, n, m)
	return pairbubble

#def Pair_Bubble(qx, qy, w):
#	return integrate.nquad(Argument, [[-inf, +inf],[-inf, +inf]])[0]

#def Dieletric_function(qx, qy, w):
#	return 1 + e**2/(2 * sqrt(qx**2 + qy**2)) * Pair_Bubble(qx, qy, w)

Qx = [ 1e-100 + 0.2 * i/50 for i in range(50)]
Qy = [ 1e-100 + 0.2 * i/50 for i in range(50)]
Kx = [-2 + 1e-100 + 4 * i/500 for i in range(500)]
Ky = [-2 + 1e-100 + 4 * i/500 for i in range(500)]

file1 = open('dieletric.dat','w')
#for kx in Qx:
#	for ky in Qy:
#		n = 1
#		Energy = E(kx, ky, n) - E(0, 0, n)		
#		fermi = f(kx, ky, 100, n, 0)
#		out = '%e    %e\n' % (Energy, fermi)
#		file1.write(out)
#		print(kx)
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
for qy in Qy:
	qx = 0
	print(qy)
	w = 1e-100
	T = 500
	Ef = 0
	n = 1
	m = 1
	q =  sqrt(qx**2 + qy**2)			
	def Argument(kx, ky):
		return Argument_real(kx, ky, qx, qy, w, T, Ef, n, m)
	pair_bubble = integrate.nquad(Argument, [[-inf, +inf],[-inf, +inf]])[0]
	dieletric_function = 2.5 + 2 * pi * pair_bubble / q		
	out = '%e    %e\n' % (q, dieletric_function)
	file1.write(out)
file1.close()

		
