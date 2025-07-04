
import numpy as np
import matplotlib.pyplot as plt
# Constants
R = 8.314  # Gas constant, J/(mol*K)
F = 96485  # Faraday constant, C/mol
T = 298  # Temperature in Kelvin
f=F/(R*T)
v = 0.01  # Potential scan rate, V/s
E_i = 0  # Initial potential, V
E_0 = 0.3  # Formal electrode potential, V
n_max=300
lambda_val=n_max/2
D= 1e-10  # in m^2/s Diffusion coefficient of oxidized species or expulsed proton
pH=0
# Create a list of E (n) values for each given n value
def E_func(var):
    E_list=[]
    for t in var:
       if t <= lambda_val:
           E=E_i - E_0+ v*t
       else:
           E=E_i-E_0 +2*v*lambda_val - v*t      
       E_list.append(E)
    return E_list

################# linear forward and backward potential scans #################

delta = 0.5
lambda_val = n_max*delta/2

def E_t(t):  
       if t <= lambda_val:
           return E_i - E_0+ 0.0591*pH+ v*t
       else:
           return E_i-E_0+ 0.0591*pH +2*v*lambda_val - v*t

################# Calculate the quasi-reversible g(n) #################

def g_quasi(n_var):
    # Initialize variables
    g= np.zeros(n_var)  # Array to store g values for each n_var

    # Loop over each n to compute g(n)
    for n in range(1, n_var):
        # Compute the summation term
        summ = 0
        for i in range(1, n-1):
            summ += np.sqrt(n - i) * (g[i + 1] - g[i])
        # Compute the right-hand side term
        param=np.sqrt(np.pi*D/(4*delta))
        k_0_term=2*k_0*np.cosh(.5*f*E_t(delta*n))/(2*k_0*np.cosh(.5*f*E_t(delta*n))+param)
        E_rev_term= param/(1+np.exp(-f*E_t(delta*n)))
        g_previous=g[n-1]-g[1]*np.sqrt(n)-summ
    
        g[n] = k_0_term*(E_rev_term + g_previous)
    return g

n_vals=np.arange(1,n_max+1)

##################### plot the quasi-reversible CV #################################

k_0_arr=[0.002,0.000002,0.0000002] # rate constants
colors=['blue','red','green']
for idx, k_0 in enumerate(k_0_arr):
        plt.plot(E_func(n_vals*delta),g_quasi(n_max)*1/.000001,color=colors[idx], linewidth=4)

plt.xticks(np.arange(-0.3, 0.51, 0.1))
plt.yticks(np.arange(-3, 4.1, 1))

plt.tick_params(axis='x', labelsize=14)  # Adjust x-axis tick label size
plt.tick_params(axis='y', labelsize=14)  # Adjust y-axis tick label size
plt.xlabel("Potential (V)", fontsize=18, labelpad=15)
plt.ylabel("Current (µA)", fontsize=18)
#plt.savefig('Fig 1.TIFF', dpi=600)