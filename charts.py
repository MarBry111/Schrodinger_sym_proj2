#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

import numpy as np
import matplotlib.pyplot as plt
import sys

kappa = 0
omega = 0
dt = 0.0001
n = 1

try:
	n = int(sys.argv[2])
	kappa = int(sys.argv[3])
	omega = int(sys.argv[4])
	print(omega)
	print("data ok")
except:
	pass

folder = "kappa"+str(kappa)+"dt"+str(dt)

rho1 = np.loadtxt(folder+"/rho_n"+str(n)+'i1k'+str(kappa)+'o'+str(omega)+".dat")
rho2 = np.loadtxt(folder+"/rho_n"+str(n)+'i2000k'+str(kappa)+'o'+str(omega)+".dat")
rho3 = np.loadtxt(folder+"/rho_n"+str(n)+'i3998k'+str(kappa)+'o'+str(omega)+".dat")

data = []

with open(folder+"/data_n"+str(n)+'k'+str(kappa)+'o'+str(omega)+".dat", "r") as f:
	for line in f:
		if line != '\n':
			data.append(line.strip("\n").split(" "))
			

data = np.array(data)

fig, axes = plt.subplots(2, 2, figsize=(18,16))
fig.suptitle('n='+str(n), fontsize=16)

axes[0, 0].plot(rho1)
axes[0, 0].plot(rho2)
axes[0, 0].plot(rho3)
axes[0, 0].set_title('rho')

axes[0, 1].plot(data[:,0], data[:,1])
axes[0, 1].set_title('Ne')
axes[1, 0].plot(data[:,0], data[:,2])
axes[1, 0].set_title('x_sr')
axes[1, 1].plot(data[:,0], data[:,3])
axes[1, 1].set_title('epsilon')

plt.savefig(folder+"/n"+str(n)+'k'+str(kappa)+'o'+str(omega)+".png")
print("save done")


