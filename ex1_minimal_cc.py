import numpy as np
import matplotlib.pyplot as plt
from gp6.gp6 import *
# import the CC data points
filename = 'https://gitlab.com/mmoresco/CCcovariance/-/raw/master/data/HzTable_MM_BC03.dat'
z_cc, Hz_cc, errHz_cc = np.genfromtxt(filename, comments = '#', usecols = (0,1,2), \
                                      unpack = True, delimiter = ',')
ref = np.genfromtxt(filename, comments = '#', usecols = (3), unpack = True, \
                    dtype = str, delimiter = ',')
gp6_cc = GP('chy', l = 1, A = 100) # random hyperparameters, choice don't matter much after optimization
gp6_cc.optimize(z_cc, Hz_cc, np.diag(errHz_cc**2)) # sets hyperparameters to optimal values given CC data
# reconstruct Hubble function and its first derivative dH/dz at Zs points
Zs = np.linspace(1e-5, 2.5, 100)
Hz_rec = gp6_cc.predict(z_cc, Hz_cc, np.diag(errHz_cc**2), Zs) # always with reference to data used for optimization
dHdz_rec = gp6_cc.predict_d1F(z_cc, Hz_cc, np.diag(errHz_cc**2), Zs)
# reconstructed H(z) plots
fig, ax = plt.subplots()
ax.errorbar(z_cc, Hz_cc, yerr = errHz_cc, fmt = 'o', color = 'red', label = 'CC data')
ax.plot(Hz_rec['z'], Hz_rec['Y'], 'b-')
ax.fill_between(Hz_rec['z'], \
                Hz_rec['Y'] - np.sqrt(Hz_rec['varY']), \
                Hz_rec['Y'] + np.sqrt(Hz_rec['varY']), \
                alpha = 0.2, color = 'blue', edgecolor = 'blue', \
                hatch = 'x', label = r'GP ($1\sigma$)', rasterized = True)
ax.set_xlim(min(Zs), max(Zs))
ax.set_xlabel('Redshift, $z$')
ax.set_ylabel('Expansion rate, $H(z)$')
ax.legend(loc = 'upper left')
fig.savefig('Hz_CC_bygp6.pdf', bbox_inches = 'tight')
# first derivative of H(z) plots
fig, ax = plt.subplots()
ax.plot(dHdz_rec['z'], dHdz_rec['Y'], 'b-')
ax.fill_between(dHdz_rec['z'], \
                dHdz_rec['Y'] - np.sqrt(dHdz_rec['varY']), \
                dHdz_rec['Y'] + np.sqrt(dHdz_rec['varY']), \
                alpha = 0.2, color = 'blue', edgecolor = 'blue', \
                hatch = 'x', label = r'GP ($1\sigma$)', rasterized = True)
ax.set_xlim(min(Zs), max(Zs))
ax.set_xlabel('Redshift, $z$')
ax.set_ylabel('Expansion rate of change, $H^\prime(z)$')
ax.legend(loc = 'lower left')
fig.savefig('dHdz_CC_bygp6.pdf', bbox_inches = 'tight')