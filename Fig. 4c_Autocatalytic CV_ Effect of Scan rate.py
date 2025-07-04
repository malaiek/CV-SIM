#Simulation of CV of Ni(OH)2 thin film electrodes at two scan rates of 2 and 5 mV/s
#Equations are described in the reference Phys. Chem. Chem. Phys., 2023, 25, 30606

import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
plt.figure(figsize=(6,6),dpi=300)

#with open('Scan-effect.csv','r') as i:
df=pd.read_csv('Scan-effect.csv')

key=np.array(df[1:],dtype=float)
xdata2=key[:,0]
ydata2=key[:,1]
plt.plot(xdata2,ydata2/0.00003142,'r', lw='2.5', label='Ni(OH)$_2$')

xdata5=key[:,2]
ydata5=key[:,3]
plt.plot(xdata5,ydata5/0.00003142,'r',lw='2.5',)

v=0.002

# Simulation with 1:3 stochiometry
Ei=0.349
x=np.linspace(1e-30,0.92,1000000)
def i3(x):
    return q*v*k3*x*(1-x)**3
def E3(x):
    return (1/k3)*(np.log(x/(1-x))+(1/(2*(1-x)**2))+(1/(1-x))-np.log(x0/(1-x0))-1/(2*(1-x0)**2)-1/(1-x0))+Ei

#FORWARD

q=0.00252
k3=770

x0=0.345
plt.plot(E3(x),i3(x)/0.00003142,'--', color='blue', lw="2.5", label='Autocatalytic Model')

#BACKWARD
q=0.002
k3=310
x0=0.81
plt.plot(E3(x),-i3(x)/0.00003142,'--', color='blue',lw='2.5')
#---------------------------------------------------------------------
v=0.005
#FORWARD
q=0.00282
k3=450

x0=0.345
plt.plot(E3(x),i3(x)/0.00003142,'--', color='blue', lw='2.5')

#BACKWARD
q=0.00212
k3=310
x0=0.81
plt.plot(E3(x),-i3(x)/0.00003142,'--', color='blue',lw='2.5')
#--------------------------------------------------------------------


plt.xlim(0.2,0.45)
plt.xticks(np.arange(min(x)+0.2, max(x)-0.45, 0.05), fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(-15,30)
plt.xlabel('Potential (V)', fontsize=18, labelpad=15)
plt.ylabel('Current (mA)',fontsize=18)
plt.legend(loc='upper right', fontsize=15)
#plt.savefig('Fig.3b.tiff', dpi=600)
