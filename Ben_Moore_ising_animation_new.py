# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 11:32:02 2020

@author: benji
"""

import numpy as np, matplotlib.pyplot as plt, matplotlib.animation as animation
%matplotlib notebook
"""
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
                if deltaE < 0:
                    s *= -1
                elif np.random.uniform(0,1) < np.exp(-deltaE*B):
                    s *= -1
                spin[x, y] = s
    return spin
"""
def E_energy(x,y): #Defining the energy due to the electric field
    s=spin[x,y]
    #enmag=s*Bfield
    spinsum = spin[(x+1)%N,y] + spin[x,(y+1)%N] + spin[(x-1)%N,y] + spin[x,(y-1)%N]
    spinprod=spinsum*s
    return -J * spinprod
    

def metro_algorithm(*args): #Applying the metropolis algorithm in a adapted way

    x = np.random.randint(N)
    y = np.random.randint(N)
    i = E_energy(x,y)
    spin[x,y] *= -1
    f = E_energy(x,y)
    deltaH = f - i
    if(np.random.uniform(0,1) > np.exp(-deltaH/T)):
        spin[x,y] *= -1

    mesh.set_array(spin.ravel())
    return mesh,


def initialise(N):  #Initialising the grid as before

    state = np.random.randint(2, size=(N,N))
    state[state==0]=-1
    return state

if __name__=="__main__": #special variables are from source file, must be changed to use as module

     N = 10 # N x N grid
     J = 1 # Exchange energy as before from energy equation
     T = 1.0 #Temperature
     minimiseE_steps = 1000 #Number of sweeps
     #opt = ''h 

     spin = initialise(N) #init_spin_config(opt) #Initial spin configuration

     #Simulation Vizualization
     fig = plt.figure(figsize=(10, 10), dpi=80)
     fig.suptitle("T = %0.1f" % T, fontsize=50)
     X, Y = np.meshgrid(range(N), range(N))
     mesh = plt.pcolormesh(X, Y, spin, cmap = plt.cm.RdBu,vmin=-1,vmax=1)
     a = animation.FuncAnimation(fig, metro_algorithm, frames = minimiseE_steps, interval = 5, blit = True)
     plt.xlim(0,N-1)
     plt.ylim(0,N-1)
     plt.show()