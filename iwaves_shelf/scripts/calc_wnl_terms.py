"""
Calculate KdV terms using WNL
"""

import numpy as np
from iwaves import IWaveModes
import matplotlib.pyplot as plt

def double_tanh(beta, z):
    return beta[0] - beta[1]*(np.tanh((z+beta[2])/beta[3])
        + np.tanh((z+beta[4])/beta[5]))



# Density parameters
#betas = [1023.7, 1.12, 105, 52, 155, 43] # ~April 5
betas = [1023.5, 1.22, 67, 55, 157, 52] # ~March 1

Z = np.arange(-250,5,5.)

rho = double_tanh(betas, Z)

iw = IWaveModes(rho,Z)

omega = 2*np.pi / (12.42*3600)
a0 = 20.

_,c1,_,_ = iw(-250,5,0)
_, Ls1 = iw.calc_steep_scales(omega,a0)

_,c2,_,_ = iw(-250,5,1)
_, Ls2 = iw.calc_steep_scales(omega,a0)
print('Ls_1, Ls_2:')
print(Ls1, Ls2)



