#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 18:17:00 2020

@author: rebeccacorley
"""

import numpy as np
import matplotlib.pyplot as plt
import emcee

# read in the table
data = np.genfromtxt('SCPUnion2.1_mu_vs_z.txt')
# pull out the redshift, magnitudes (mm) and magnitude estimated errors (dm)
redshift = data.T[1]
mm = data.T[2]
dm = data.T[3]
dlabel = 'SCP_2.1'

# mm is the "distance modulus" 
# from which we can compute the distance in parsecs (pc)
dpc = 10.**(mm/5.+1.)
# and in megaparsecs (Mpc)
dMpc = dpc / 10.**6
# and the error on that distance:
dMe = 10.**((mm+dm)/5.+1.-6.) - dMpc

# Likelihood statement for model similar to FRW with OmM = 0.43


def model(c_H0, redshift):
    return c_H0*(redshift)*np.sqrt(1+redshift)

def log_likelihood(parameters):
    c_H0 = parameters
    exp_model = model(c_H0, redshift)
    sigma = dMe
   
    X2 = ((dMpc - exp_model)/sigma)**2
    L_i = np.sum(-0.5*np.log((2*np.pi*(sigma)**2)) -  X2/2.0)
   
    return L_i


nwalkers, ndim = 500,1
nsteps = 1000

#randomly seed intial distribution

c_H0_0 = np.random.uniform(67, 74, nwalkers)

#dMpc_0 = np.random.uniform(0,30, nwalkers)

#initial_state = np.column_stack((c_H0_0, dMpc_0))
initial_state = c_H0_0.reshape(nwalkers, 1)
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood)
sampler.run_mcmc(initial_state, nsteps)

walker_chains = sampler.flatchain
np.savetxt('hubble_chain_2.txt', walker_chains)

chain = np.genfromtxt('hubble_chain_2.txt')
chain = chain.reshape(1,1000,500)

for c in chain:
    plt.plot(c)

plt.tight_layout()
plt.savefig('SCP_FRW_chains.pdf')
plt.xlabel('Step Number')
plt.show()

H0 = 3e5/4500
sH0 = 'H0 = '+str(round(H0))+' km/s/Mpc'

plt.plot(redshift, dMpc, '.')
plt.plot(redshift, (4500*redshift*np.sqrt(1+redshift)),label="d = c/H0*z*sqrt(1+z)")
plt.xlabel(' z')
plt.ylabel('Luminosity distance (Mpc)')
plt.text(0.1,9000,sH0,fontsize = 12)
plt.grid(b=True,which='both')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('SCP_FRW_graph.pdf')
plt.show()


print(H0)