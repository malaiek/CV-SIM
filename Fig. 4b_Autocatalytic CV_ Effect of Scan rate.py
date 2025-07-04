#Simulation of the first and second CV cycles of Ni(OH)2 thin film electrodes
#Equations are described in the reference Phys. Chem. Chem. Phys., 2023, 25, 30606

import numpy as np
import matplotlib.pyplot as plt
import csv
plt.figure(figsize=(5,5),dpi=300)



v=0.002


with open('Ni-2mv.csv','r') as i:
    xydata = list(csv.reader(i,delimiter=","))
key=np.array(xydata[1:],dtype=float)
xdata=key[:,0]
ydata=key[:,1]
plt.plot(xdata,ydata/0.00003142,'r', label='Experiment', lw='3')

# Simulation with 1:3 stochiometry
Ei=0.345
x=np.linspace(1e-30,0.92,1000000)
def i3(x):
    return q*v*k3*x*(1-x)**3
def E3(x):
    return (1/k3)*(np.log(x/(1-x))+(1/(2*(1-x)**2))+(1/(1-x))-np.log(x0/(1-x0))-1/(2*(1-x0)**2)-1/(1-x0))+Ei
q=0.0042
#FORWARD
k3=800
x0=0.00000005
plt.plot(E3(x),i3(x)/0.00003142,'--', color='b', lw='3', label='Autocatalytic Model')

#BACKWARD
q=0.0021
k3=320
x0=0.81
plt.plot(E3(x),-i3(x)/0.00003142,'--', color='b',lw='3')

#FORWARD
q=0.0026
k3=800

x0=0.367
plt.plot(E3(x),i3(x)/0.00003142,'--', color='blue', lw='3')

#BACKWARD
q=0.0026
k3=320
x0=0.81
plt.plot(E3(x),-i3(x)/0.00003142,'--', color='blue',lw='3')


plt.xticks(np.arange(min(x)+0.2, max(x)-0.45, 0.05))
plt.ylim(-10,30)
plt.xlim(0.2,0.45)
plt.xlabel('Potential (V vs. Ag/AgCl)',fontsize=12)
plt.ylabel('Current (mA/$\mathregular{cm^2}$)',fontsize=12)
plt.legend(loc='upper left')
plt.legend()
#plt.savefig('Fig.3a.tiff', dpi=600)


