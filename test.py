from math import *

EF = 3 # Fermi energy in eV
hbar = 0.65821e-15#Dirac's constant in eV . s
h = 4.13566743e-15#Planck's constant in eV.s
e = 1.60217662e-19 #Quantum of charge in C
mx = (hbar**2)/(2 * 3.775 ) 
my = (hbar**2)/(2 * 19.564 ) 

def plus_f(q, w):
	if (w/q + q/2)**2 - 2*EF*my/hbar**2 != 0 and 0 < (2*EF*mx/hbar**2 - (w/q + q/2)**2)*((w/q + q/2)**2 - 2*EF*my/hbar**2):
		value = sqrt((2*EF*mx/hbar**2 - (w/q + q/2)**2)/((w/q + q/2)**2 - 2*EF*my/hbar**2))
		return value
	else:
		return 0

def minus_f(q, w):
	if (w/q - q/2)**2 - 2*EF*my/hbar**2 != 0 and 0 < (2*EF*mx/hbar**2 - (w/q - q/2)**2)*((w/q - q/2)**2 - 2*EF*my/hbar**2):
		value = sqrt((2*EF*mx/hbar**2 - (w/q - q/2)**2)/((w/q - q/2)**2 - 2*EF*my/hbar**2 ))
		return value
	else:
		return 0

frequencies = [1e-6 + 2 * i/1000 for i in range(1000)]
momenta = [1e-6 + 2 * i/1000 for i in range(1000)]
#qx = 0.5
#qy = 0.5
#q = sqrt((qx**2 )*(hbar**2)/mx + (qy**2 )*(hbar**2)/my)

file1 = open('test.dat','w')
for w in frequencies:
	q = 1
	positive = plus_f(q, w)
	negative = minus_f(q, w)
	tg = tan(w)	
	out = '%e    %e    %e    %e\n' % (w, positive, negative, tg)
	file1.write(out)
file1.close()
