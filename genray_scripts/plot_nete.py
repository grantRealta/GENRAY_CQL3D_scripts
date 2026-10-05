###
# Plots the electron and ion densities and temperatures according to the genray input file
# The plot_profs.py cql3d script is superior, but this can be used as a sanity check
###

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import f90nml
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__);dname = os.path.dirname(abspath);os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__));parentdir = os.path.dirname(currentdir);sys.path.append(parentdir)

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 17)
plt.rc('figure', titlesize = 16)

genrayPath = f'{targetDir}/genray.in'
nml = f90nml.read(genrayPath)

temtab = nml['temtab']['prof']
dentab = nml['dentab']['prof']

numSpecies = len(nml['species']['charge'])


T_e = np.zeros(int(len(temtab)/numSpecies))
#T_i = np.copy(T_e)
n_e = np.copy(T_e)
#n_i = np.copy(T_e)

#read in the profiles
for i in range(0,len(T_e)):
    T_e[i] = temtab[numSpecies*i]
    #T_i[i] = temtab[numSpecies*i+1]
    n_e[i] = dentab[numSpecies*i]
    #n_i[i] = dentab[numSpecies*i+1]

tescal = 1
try:
    tescale = nml['plasma']['temp_scale(1)']
except:
    pass
nescal = 1
try:
    nescale = nml['plasma']['den_scale(1)']
except:
    pass

#apply any scalings
T_e = T_e*tescal
#T_i = T_i*tiscal
n_e = n_e*nescal
#n_i = n_i*niscal
rhosHelper = np.arange(1,len(T_e)+1,1)
rhos = (rhosHelper-1)/(len(T_e)-1)

fig, ax = plt.subplots(figsize = (8*.8,6*.8))
ax.plot(rhos, np.array(n_e/1e19), label = r'$n_e$', linestyle = 'dashed', color = 'k', linewidth = 3)
#ax.plot(rhos, T_i, label = r'$T_i$', linewidth = 3, color = 'b')
ax2 = ax.twinx()
ax2.plot(rhos, T_e*1e3, label = r'$T_e$', linewidth = 3, color = 'g')
#ax.plot(rhos, np.array(n_i/1e19), label = r'$n_i$', linestyle = 'dotted', color = 'k', linewidth = 3)

#ax.legend(loc = 'best', fontsize = 20, ncol=2)

ax.set_xlabel(r'$\rho_{pol}$', fontsize = 20)
ax.set_xlim([0,1])
#ax.set_ylim([0,8])
ax.set_ylabel(r'n$_j$ (10$^{19}$/m$^{-3}$)')
ax2.set_ylabel(r'T$_j$ (eV)')
ax.grid()
fig.tight_layout()
#fig.tight_layout(rect=[0, 0.0, 1, 1.05])
plt.show()

