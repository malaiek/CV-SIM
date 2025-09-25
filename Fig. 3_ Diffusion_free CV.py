###ploting CV curves for the diffusion-free EPT
### Underlying model development is described in reference ChemElectroChem 2023, 10, e202201118

import numpy as np
import matplotlib.pyplot as plt
R, T, F, n, b, q, v, E_0 = 8.314, 298, 96485.3, 1, 1, 1, 0.01,.15
r=R*T/F
f=38.9
fig,axs=plt.subplots(1,2, figsize=(12,6))

### plotting Fig. a
def E_rxn_cntr(x,g,E_0,a_ion): # Nernst equation modfied with ineraction parameter g
    return E_0+(b/n)*r*np.log(a_ion)+r/n*np.log((1-x)/x)+r/n*(g/2)*(2*x-1)
x_data=np.linspace(0.01,0.99,1000)

def i_rxn_cntr(x,g):
    return f*q*v*(x*(1-x)/(g*x*(1-x)-1))
axs[0].plot(E_rxn_cntr(x_data,0,.15,1),i_rxn_cntr(x_data,0),'b',linewidth='4', label='Homogeneous O/R')
axs[0].plot(E_rxn_cntr(x_data,0,.15,1),-i_rxn_cntr(x_data,0),'b',linewidth='4')
axs[0].plot(E_rxn_cntr(x_data,2,.15,1),i_rxn_cntr(x_data,2),'r',linewidth='4',label='Heterogeneous O/R')
axs[0].plot(E_rxn_cntr(x_data,2,.15,1),-i_rxn_cntr(x_data,2),'r',linewidth='4')

### plotting Fig. b

x_data=np.linspace(0.0001,0.9999,1000)

#forward reaction
def E_hyster_ox(x,gamma,alpha): # alpha and gamma are parameters that account for reaction path
    return E_0+r*(np.log((1-x)/x)-(gamma/4)*(1-2*x)-gamma*alpha*np.exp(alpha*x))
def i_hyster_ox(x,gamma,alpha):
    return -(f*q*v)*x*(1-x)/((gamma*(.5-alpha**2*np.exp(alpha*x)))*x*(1-x)-1)

#backward reaction
def E_hyster_red(x,gamma,alpha):
    return E_0+r*(np.log((1-x)/x)-(gamma/4)*(1-2*x)-gamma*(.5-alpha)*np.exp((.5-alpha)*x))
def i_hyster_red(x,gamma,alpha):
    return (f*q*v)*x*(1-x)/((gamma*(.5-(.5-alpha)**2*np.exp((.5-alpha)*x)))*x*(1-x)-1)

axs[1].plot(E_hyster_ox(x_data,2,.01),i_hyster_ox(x_data,2,.01),'b',linewidth='4',label='Hetero O/R: weak hysteresis')
axs[1].plot(E_hyster_red(x_data,2,.01),i_hyster_red(x_data,2,.01),'b',linewidth='4')
axs[1].plot(E_hyster_ox(x_data,6,.01),i_hyster_ox(x_data,6,.01),'r',linewidth='4',label='Hetero O/R: strong hysteresis')
axs[1].plot(E_hyster_red(x_data,6,.01),i_hyster_red(x_data,6,.01),'r',linewidth='4')


#### fonts and labels####

axs[0].set_xlabel("Potential (V)", fontsize=18, labelpad=15)
axs[0].set_ylabel("Current (mA)", fontsize=18)
axs[0].set_xticks(np.arange(0,0.3,0.05))
axs[0].set_yticks(np.arange(-0.25,.35,0.1))
axs[0].tick_params(axis='x', labelsize=15)  
axs[0].tick_params(axis='y', labelsize=15)

axs[1].set_xlabel('Potential(V)', fontsize=18, labelpad=15)
axs[1].set_ylabel('Current (mA)', fontsize=18)
axs[1].tick_params(axis='x', labelsize=15)  
axs[1].tick_params(axis='y', labelsize=15)
axs[1].set_xticks(np.arange(-0.2,0.45,0.1))
axs[1].set_yticks(np.arange(-0.2,0.5,0.1))

axs[0].legend(fontsize=12)
axs[1].legend(loc='best',fontsize=10)
plt.tight_layout()
#plt.savefig('Fig.3.tiff',dpi=600)

