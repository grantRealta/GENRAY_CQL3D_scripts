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
import helperFunctions as helper

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 18)
plt.rc('figure', titlesize = 20)
plt.rc('legend', fontsize = 14)
plt.rcParams['lines.linewidth'] = 2

genrayPath = f'{targetDir}/genray.in'
nml = f90nml.read(genrayPath)

dentab = nml['dentab']['prof']

freq = 35e9#nml['wave']['frqncy']
omega = freq*2*np.pi
numSpecies = len(nml['species']['charge'])

n_e = np.zeros(int(len(dentab)/numSpecies))

for i in range(0,len(n_e)):
    n_e[i] = dentab[numSpecies*i]

nescal = 1
try:
    nescale = nml['plasma']['den_scale(1)']
except:
    pass

n_e = n_e*nescal
print(f'max of ne: {np.nanmax(n_e)}')
rhosHelper = np.arange(1,len(n_e)+1,1)
rhos = (rhosHelper-1)/(len(n_e)-1)


w_pe2 = helper.get_wpe2(n_e)
w_pD2 = helper.get_wpD2(n_e)

Nparas = [0,.5,1,1.5]
B = 1.25

fig,ax = plt.subplots()

for i, Npara in enumerate(Nparas):
    W_ce = helper.get_Wce(B)
    W_cD = helper.get_WcD(B)

    S = 1 - (w_pe2/(omega**2 - W_ce**2)) - (w_pD2/(omega**2 - W_cD**2))
    D = (W_ce/omega)*(w_pe2/(omega**2-W_ce**2)) + (W_cD/omega)*(w_pD2/(omega**2-W_cD**2))
    P = 1 - (w_pe2/omega**2) - (w_pD2/omega**2)

    C_0 = P*((Npara**2-S)**2-D**2)
    C_2 = (Npara**2 - S)*(P+S) + D**2
    C_4 = S

    Nperp2_minus = (-C_2 - np.sqrt(C_2**2-4*C_4*C_0))/(2*C_4)
    Nperp2_plus = (-C_2 + np.sqrt(C_2**2-4*C_4*C_0))/(2*C_4)

    ax.plot(n_e, Nperp2_plus, label = f'$N_{{||}}$ = {Npara}')
        #ax.plot(rhos, Nperp2_minus, label = f'-root, N|| = {Npara}, B = {B} T', linestyle = 'dotted')

ax.axhline(0, color = 'k', linestyle = 'dashed')

ax.set_title(f'f = {freq/1e9:.0f} GHz, B = {B} T')
ax.legend(loc = 'best')
ax.set_ylabel(r'$N_{\perp}^2$')
ax.set_xlabel(r'$n_e$ ($m^{-3}$)')
fig.tight_layout()
plt.show()