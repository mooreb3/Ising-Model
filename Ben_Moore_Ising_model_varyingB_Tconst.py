# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 10:57:37 2020

@author: benji
"""


#constant temperature and changing magnetic field
#magnetisation and susceptibility do not change with B field once equilib is reached
import numpy as np
import matplotlib.pyplot as plt

print('Benjamin Moore 17336557') 
print('Constant Temperature and Magnetic Field Varying')

def initialise(N):  

    state = np.random.randint(2, size=(N,N))
    state[state==0]=-1
    return state
#Initialising lattice to have random spin (plus or minus 1)
def metro_algorithm(spin, B, Bfield): #metropolis algorithm as outlined in method
    for i in range(N):
        for j in range(N):
                x = np.random.randint(N)
                y = np.random.randint(N)
                s=spin[x,y]
                enmag=s*Bfield
                spinsum = spin[(x+1)%N,y] + spin[x,(y+1)%N] + spin[(x-1)%N,y] + spin[x,(y-1)%N]
                spinprod=spinsum*s
                E_energy=J*spinprod
                total_E=E_energy+enmag #suming electrical and magnetic contributions, magnetic is zero
                deltaE = 2*total_E
                if deltaE < 0:  #First criteria to decide whether to flip
                    s *= -1
                elif np.random.uniform(0,1) < np.exp(-deltaE*B): #2nd criterial; Flip due to probability
                    s *= -1
                spin[x, y] = s
    return spin


def Energy_of_sys(spin,Bfield):         #Energy of system defined as outlined in theory, with the electronic energy being divided by 4 to account for the neighbours(divided by number of sites later)
    energy = 0
    for i in range(len(spin)):
        for j in range(len(spin)):
            s = spin[i,j]
            enmag=s*Bfield
            spinsum = spin[(i+1)%N, j] + spin[i,(j+1)%N] + spin[(i-1)%N, j] + spin[i,(j-1)%N]  #modulo funtion allows the shape to be periodic by linking the edge inputs
            energy += -spinsum*J*s/4-enmag
    return energy
 
def Magnetisation(spin,dB):
    mag = np.sum(spin)
    return mag
J=1
N       = 5         #  size of the lattice, N x N
minimiseE_steps = 1000      #  number of MC sweeps for equilibration
ave_sweeps = 1000       #  number of MC sweeps for calculation
no_dB=30
Bfield = np.linspace(0.0,10.0,no_dB)  
T=1
E,M,C,X = np.zeros(no_dB), np.zeros(no_dB), np.zeros(no_dB), np.zeros(no_dB)
normalise_1  = 1.0/(ave_sweeps*N*N)
normalise_2 = 1.0/(ave_sweeps*ave_sweeps*N*N)
#Normalisation factor for average energy and magnetisation as outlined before


for a in range(no_dB):
    E1 = M1 = E2 = M2 = 0
    spin = initialise(N)
    B=1.0/T
    B2=B*B
    dB=1.0*Bfield[a]

 
    for i in range(ave_sweeps):
        metro_algorithm(spin, B, dB)          
        Ene = Energy_of_sys(spin, dB)     # calculate the energy
        Mag = Magnetisation(spin, dB)        # calculate the magnetisation
 
    E1 = E1 + Ene
    M1 = M1 + Mag
    M2 = M2 + Mag*Mag
    E2 = E2 + Ene*Ene
 
    E[a] = normalise_1*E1
    M[a] = normalise_1*M1
    C[a] = (normalise_1*E2 - normalise_2*E1*E1)*B2
    X[a] = (normalise_1*M2 - normalise_2*M1*M1)*B
    

 
#f = plt.figure(figsize=(18, 10)); # plot the calculated values   
 
 
plt.figure()
plt.plot(Bfield, E, marker='o', color='IndianRed')
plt.xlabel("B field (T)", fontsize=20);
plt.ylabel("Energy(J) ", fontsize=20);         plt.axis('tight');
 
 
plt.figure()
plt.plot(Bfield, abs(M), marker='o', color='RoyalBlue')
plt.xlabel("B field(T)", fontsize=20);
plt.ylabel("Magnetization (A/m) ", fontsize=20);   plt.axis('tight');
 
plt.figure()
plt.plot(Bfield, C, marker='o', color='IndianRed')
plt.xlabel("B field(T)", fontsize=20); 
plt.ylabel("Specific Heat (J/K) ", fontsize=20);   plt.axis('tight');  
 
plt.figure()
plt.plot(Bfield, X, marker='o', color='RoyalBlue')
plt.xlabel("B field (T)", fontsize=20);
plt.ylabel("Susceptibility", fontsize=20);   plt.axis('tight');


