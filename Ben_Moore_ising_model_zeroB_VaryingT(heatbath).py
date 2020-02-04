# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 16:44:43 2020

@author: benji
"""
#Ising model for zero magnetic field and changing temperature

import numpy as np
import matplotlib.pyplot as plt

print('Benjamin Moore 17327505')
print('Varying Temperature at B=0')  
def initialise(opt):  

    state = np.random.randint(2, size=(N,N))
    state[state==0]=-1
    return state
#Initialising the random lattice with a random grid of numbers 1 and -1 representing spin
enmag=0
def metro_algorithm(spin, B):
    probvals = [np.exp(-abs(enmag)*2*B),np.exp(-(2+abs(enmag))*2*B),np.exp(-(4+abs(enmag))*2*B), np.exp(-(6+abs(enmag)*2*B))]
    #only 3 values available for probability-saves computational time
    #n=1
    for i in range(N):
        for j in range(N):
                x = np.random.randint(N)
                y = np.random.randint(N)
                s=spin[x,y]
                #enmag=s*magf
                spinsum = spin[(x+1)%N,y] + spin[x,(y+1)%N] + spin[(x-1)%N,y] + spin[x,(y-1)%N] #modulo funtion allows the shape to be periodic by linking the edge inputs
                spinprod=spinsum*s
                E_energy=J*spinprod
                #total_E=E_energy-enmag suming electrical and magnetic contributions, magnetic is zero
                probi=int((E_energy/2-enmag/2))
                if E_energy < 0:
                    s *= -1
                elif np.random.uniform(0,1) < probvals[probi]:
                    s *= -1
                spin[x, y] = s
    return spin

def Energy_of_sys(spin):    #Energy of system defined as outlined in theory, with the electronic energy being divided by 4 to account for the neighbours
    energy = 0
    for i in range(len(spin)):
        for j in range(len(spin)):
            s = spin[i,j]
            enmag=s*magf
            spinsum = spin[(i+1)%N, j] + spin[i,(j+1)%N] + spin[(i-1)%N, j] + spin[i,(j-1)%N]
            energy += -spinsum*J*s/4-enmag
    return energy
 
 
def Magnetisation(spin):
    '''Magnetization of a given configuration'''
    mag = np.sum(spin) #Defining function of total spins, will be divided by number of sites later
    return mag


magf=0.0
J=1
tempstep      = 100        #steps in the temperature range
N       = 5       #NxN grid
minimiseE_steps = 1000     #number of sweeps to reach minimum energy
ave_sweeps = 1000   #Monte Carlo sweeps used to get average of observables from final n sweeps
T = np.linspace(1.0, 5.0, tempstep); #defining the temperature range
n1, n2  = 1.0/(ave_sweeps*N*N), 1.0/(ave_sweeps*ave_sweeps*N*N)
# divide by number of samples (in eq for ave energy and magnetism) and by system size to get average values 

 #making empty observables list to store
E= np.zeros(tempstep)         
M= np.zeros(tempstep)
C= np.zeros(tempstep)
X = np.zeros(tempstep)






for tt in range(tempstep):
    E1 = M1 = E2 = M2 = 0
    spin = initialise(N)
    iT=1.0/T[tt]; iT2=iT*iT;
    
    for i in range(ave_sweeps):
        metro_algorithm(spin, iT)          
        Ene = Energy_of_sys(spin)     # calculat the energy
        Mag = Magnetisation(spin)        #and magnetisation of system per tempstep to be averaged
 
        E1 = E1 + Ene   
        M1 = M1 + Mag
        M2 = M2 + Mag*Mag
        E2 = E2 + Ene*Ene
 
    E[tt] = n1*E1 #Normalising and averaging the 4 observables
    M[tt] = n1*M1
    C[tt] = (n1*E2 - n2*E1*E1)*iT2
    X[tt] = (n1*M2 - n2*M1*M1)*iT
    

 
#Plotting
plt.figure()
plt.plot(T, E, marker='o', linestyle='None')
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Energy ", fontsize=20);         plt.axis('tight');
 
plt.figure()
plt.plot(T, abs(M), marker='o',linestyle='None', color='RoyalBlue')
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight');
 
plt.figure()
plt.plot(T, C, marker='o',linestyle='None', color='IndianRed')
plt.xlabel("Temperature (T)", fontsize=20); 
plt.ylabel("Specific Heat ", fontsize=20);   plt.axis('tight');  
 
plt.figure()
plt.plot(T, X, marker='o',linestyle='None', color='RoyalBlue')
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Susceptibility", fontsize=20);   plt.axis('tight'); 

