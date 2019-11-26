# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 05:58:18 2019

@author: Usuario
"""

import numpy as np;


kb= 1.38*10**(-23)
V=1000
N=500
f=3

K, U, E = np.loadtxt("results_energies.txt", unpack=True)

K_inv, Phi_VK_inv, Phi_V, sqrd_Phi_VK_inv, Phi_VV= np.loadtxt("mean_values.txt", unpack= True)



K_m= np.mean(K)
U_m= np.mean(U)
E_m= np.mean(E)

T= 2/(f*kb)*K_m
P= N/V*kb*T-Phi_V
Cv= kb/(1+(2/f-1)*K_m*K_inv)
alpha_E= kb/((1-2/f)*K_m*Phi_VK_inv- Phi_V)
gamma= N*kb/Cv+V*(f/2-1)*(Phi_V*K_inv- Phi_VK_inv)
ks=1/( N*kb*T/V*(1+2*gamma-N*kb/Cv)+V*Phi_VV-V*(f/2-1)*(sqrd_Phi_VK_inv-2*Phi_V*Phi_VK_inv+Phi_V**(2)*K_inv))


alpha_s= -1/(gamma*T)
kT= 1/(1/ks-T*Cv/V*gamma**(2))
Cp= Cv*(kT/ks)
alpha_P=Cv/V*gamma*kT



print("T= ", T, "\n",  "P= ", P, "\n", "Cv= ", Cv, "\n",  "alpha_E= ", alpha_E, "\n",  "gamma= ", gamma, "\n", "K_s= ",  ks)
print("alpha_S= ",  alpha_s,"\n", "K_T= ", kT, "\n",  "C_P= ", Cp, "\n",  "alpha_P= ", alpha_P)