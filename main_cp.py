# -*- coding: utf-8 -*- 

#to use when we know for which omega we should use it

import numpy as np
import sys
import os

def generate_H(psi, x_k, delta_x, kappa, omega, tau, N):
	H = np.zeros(shape=(1,N+1))
	H[:,1:N] = -1*(psi[:,:N-1]+psi[:,2:]-2*psi[:,1:N])/delta_x/delta_x/2 + kappa*(x_k[:,1:N]-1/2)*psi[:,1:N]*np.sin(omega*tau)
	return H

percent = 1.0 # wczytywane
const = 0 #to zaraz wybierzemy
kappa = 5
tau = 0

N = 100
delta_x = 1./N
x_k = np.array(range(N+1)).reshape(1,N+1)*delta_x

n = 1 #1 4 9

number_of_steps = 4000
delta_tau = 0.0001
psi_I = np.zeros(shape=(1,N+1))
psi_R = np.sqrt(2)*np.sin(np.pi*x_k*n)

try:
	percent = float(sys.argv[2])/100
	print(percent)
except:
	print("except data")

folder = "kappa"+str(kappa)#+"dt"+str(delta_tau)

omega = percent*(const/2)*np.pi*np.pi 

H_R = generate_H(psi=psi_R, x_k=x_k, delta_x=delta_x, kappa=kappa, omega=omega, tau=tau, N=N)
H_I = generate_H(psi=psi_I, x_k=x_k, delta_x=delta_x, kappa=kappa, omega=omega, tau=tau, N=N)

try:
	os.mkdir(folder)
except:
	print("except folder")

path_dat = folder+"/data_n"+str(n)+'k'+str(kappa)+'p'+str(percent)+".dat"

try: 
	f = open(path_dat, "a")
	f.close()
except:			
	with open(path_dat, "w") as f:
		f.write('')

epsilon = 0
for i in range(number_of_steps):
	tau = i*delta_tau
	psi_R = psi_R + H_I*delta_tau/2
	H_R = generate_H(psi=psi_R, x_k=x_k, delta_x=delta_x, kappa=kappa, omega=omega, tau=tau, N=N)
	psi_I = psi_I - H_R*delta_tau
	H_I = generate_H(psi=psi_I, x_k=x_k, delta_x=delta_x, kappa=kappa, omega=omega, tau=tau, N=N)
	psi_R = psi_R + H_I*delta_tau/2
	
	epsilon_new = delta_x*np.sum( psi_R*H_R + psi_I*H_I )
	if epsilon_new > epsilon: 
		epsilon = epsilon_new

with open(path_dat, "a") as f:
			f.write(str(omega)+' '+str(epsilon))
			f.write('\n')

print("finish")
