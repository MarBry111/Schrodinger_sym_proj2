# -*- coding: utf-8 -*- 

import numpy as np
import sys
import os

#kilka omeg

#nadowan czstka 1D przestrzen i jednorodne pole o zmiennej amplitudzie oscylujacej w czasie
#ewolucja f falowej elektrony i pochlanianie i oddawnaie energii elektrycznej

#hamiltionian 3, 2, 13 - suma stanow stacjonarnych
#przyblizanie - 14 z dwuch stanow stacjonarnych

#strona 4 - parametryzacja
'''
X - pd 0 do L opis przestrzeni 1D
16
x - stosunek wartosci w przestrzeni X do L (L - wymiar ukladu)  - polozenie
tau - parametryzuje stala plancka, czas wlasciwy, mase i wymiar ukladu  - czas

PARAMETRYZACJA POZBAWIA NAS JEDNOSTEK FIZYCZNYCH
19 - parametrezacja pola
kappa - 
omega -

ALGORYTM NUMERYCZNY - strona 5
- komorki od 0 do N
29 - warunki poczatkowe do rozw rownania Schrodingera
30 - hamiltionian (2 czesc - zmienne pole)
	czesc rzeczywista i urojona taka sama
r schrodingera - czesc rzeczywista i urojona
ten sam typ algorytmu co dla argonu
dzialamy na f stanu - f falowej

norma, polozenie i energia
gestosc prawdopodobienstwa - zalezne nie od czasu a polozenia - w jednostkach obliczneiach
k - numer siatki obliczeniowej

kropk calkowania delta tau
zbadac we wzorcowym rozwiazaniu - bez pola
potem z wlaczonym polem, obszar rezonansowy

tau N x E ro
N(tau) x(tau)
E(tau) ro(k) - 3 serie (na poczatku , w srodku i na koniec)

30 kappa=0 omega=9 - brak zabuurzen zewnetrznego pola
'''

#gnuplot - plot zamiast splot

def generate_H(psi, x_k, delta_x, kappa, omega, tau, N):
	H = np.zeros(shape=(1,N+1))
	H[:,1:N] = -1*(psi[:,:N-1]+psi[:,2:]-2*psi[:,1:N])/delta_x/delta_x/2 + kappa*(x_k[:,1:N]-1/2)*psi[:,1:N]*np.sin(omega*tau)
	return H

const = 4
kappa = 5
tau = 0

N = 100
delta_x = 1./N
x_k = np.array(range(N+1)).reshape(1,N+1)*delta_x

n = 1 #1 4 9

number_of_steps = 600000
delta_tau = 0.0001

psi_I = np.zeros(shape=(1,N+1))
psi_R = np.sqrt(2)*np.sin(np.pi*x_k*n)

try:
	const = int(sys.argv[2])
	#kappa = int(sys.argv[3])
	print(const)
except:
	print("except data")

folder = "kappa"+str(kappa)#+"dt"+str(delta_tau)

omega = (const/2)*np.pi*np.pi 

H_R = generate_H(psi=psi_R, x_k=x_k, delta_x=delta_x, kappa=kappa, omega=omega, tau=tau, N=N)
H_I = generate_H(psi=psi_I, x_k=x_k, delta_x=delta_x, kappa=kappa, omega=omega, tau=tau, N=N)

try:
	os.mkdir(folder)
except:
	print("except folder")

path_dat = folder+"/data_n"+str(n)+'k'+str(kappa)+'o'+str(const)+".dat"

with open(path_dat, "w") as f:
	f.write('')

#rho stuff
path_rho = "rho/rho_"+str(n)+'k'+str(kappa)+'o'+str(const)+".dat"

try: 
	f = open(path_rho, "a")
	f.close()

except:			
	f = open(path_rho, "w")
	f.close()

T = 2*np.pi/omega
print(T)
iterat = 0

for i in range(number_of_steps):
	tau = i*delta_tau
	psi_R = psi_R + H_I*delta_tau/2
	H_R = generate_H(psi=psi_R, x_k=x_k, delta_x=delta_x, kappa=kappa, omega=omega, tau=tau, N=N)
	psi_I = psi_I - H_R*delta_tau
	H_I = generate_H(psi=psi_I, x_k=x_k, delta_x=delta_x, kappa=kappa, omega=omega, tau=tau, N=N)
	psi_R = psi_R + H_I*delta_tau/2
	
	if i%300 ==0:
		with open(path_dat, "a") as f:
			Ne = delta_x*np.sum( psi_R**2 + psi_I**2 )
			x_sr = delta_x*np.sum( x_k*(psi_R**2 + psi_I**2) )
			epsilon = delta_x*np.sum( psi_R*H_R + psi_I*H_I )
			f.write(str(tau+delta_tau+i*delta_tau)+' '+str(Ne)+' '+str(x_sr)+' '+str(epsilon))
			f.write('\n\n')
		
		rho_k = psi_R[:,::2]**2 + psi_I[:,::2]**2
		with open(path_rho, "a") as f:
			for j in rho_k[0]:
				f.write(str(j))
				f.write("\n")
			f.write("\n\n")
		
		print(float(i)/number_of_steps)
		
	if (i == 1 or i == number_of_steps/2 or i ==number_of_steps-2) and iterat<3:
		rho_k = psi_R[:,::2]**2 + psi_I[:,::2]**2
		np.savetxt(folder+"/rho_n"+str(n)+'i'+str(i)+'k'+str(kappa)+'o'+str(const)+".dat", rho_k)
		iterat+=1

os.system("python charts.py -i "+str(n)+" "+str(kappa)+" "+str(const))

print("finish")
