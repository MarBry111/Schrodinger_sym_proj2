import numpy as np
import matplotlib.pyplot as plt

kappa = 5
omega = 3
tau = 0
dat_i = 0.0001

folder = "kappa"+str(kappa)+"omega"+str(omega)
n = 1

rho1 = np.loadtxt(folder+'/rho_n'+str(n)+'i1.dat')
rho2 = np.loadtxt(folder+'/rho_n'+str(n)+'i2000.dat')
rho3 = np.loadtxt(folder+'/rho_n'+str(n)+'i3998.dat')

data = []

with open(folder+"/data_n"+str(n)+".dat", "r") as f:
	for line in f:
		if line != '\n':
			data.append(line.strip("\n").split(" "))
			

data = np.array(data)

fig, axes = plt.subplots(2, 2)
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

plt.savefig(folder+"/n"+str(n)+".png")
print("save done")


