###
# Plots the electron and ion densities and temperatures according to the cql3d input file
# Either plots them all on one plots or a two axis plots (one for T, one for n)
###

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__);dname = os.path.dirname(abspath);os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__));parentdir = os.path.dirname(currentdir);sys.path.append(parentdir)

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

import getInputFileDictionary
inputFileDict = getInputFileDictionary.getInputFileDictionary('cql3d')

import numpy as np
import matplotlib.pyplot as plt
import helperFunctions as helper

splitTandN = False

if not splitTandN:
    nrows = 1
else:
    nrows = 2

rho_pol, n_e = helper.getCQLne()
rho_pol, T_e = helper.getCQLTe()
rho_pol, n_D = helper.getCQLnD()
rho_pol, T_D = helper.getCQLTD()

dischargeNumber = os.getcwd().split('/')[-1]
if len(dischargeNumber) != 6:
    dischargeNumber = os.getcwd().split('/')[-2]

#####Setup and do plotting#####
plt.rc('xtick', labelsize = 15)
plt.rc('ytick', labelsize = 15)
plt.rc('axes', labelsize = 17)
plt.rc('figure', titlesize = 16)
plt.rc('legend',fontsize=16)

if not splitTandN:
    figsize = (6.4,4.5)#(6.4,4.8)
else:
    figsize = (6.4,5.5)

fig,axes = plt.subplots(nrows = nrows, figsize=figsize)

if not splitTandN:
    axes.set_xlim([0,1])
    axes2 = axes.twinx()

    axes.plot(rho_pol, n_e/1e19, label = r'$n_e$', lw = 3, color = 'green')
    axes2.plot(rho_pol, T_e, label = r'$T_e$', lw = 3,linestyle = 'dashed', color = 'green')
    axes.plot([.2,.4],[-1,-2], label = r'$T_e$', lw = 3,linestyle = 'dashed', color = 'green')

    axes.plot(rho_pol, n_D/1e19, label = r'$n_D$', lw = 3, color = 'royalblue')
    axes2.plot(rho_pol, T_D, label = r'$T_D$', lw = 3,linestyle = 'dashed', color = 'royalblue')
    axes.plot([.2,.4],[-1,-2], label = r'$T_D$', lw = 3,linestyle = 'dashed', color = 'royalblue')

    axes.set_ylabel(r'Density ($10^{19}m^{-3}$)')
    axes2.set_ylabel(r'Temperature (keV)')
    axes.set_xlabel(r'$\rho_{pol}$')
    axes.grid()
    axes.legend(ncol = 2)
    axes.set_ylim([0,10])
    axes2.set_ylim([0,10])
    

else:

    axes[0].plot(rho_pol, n_e/1e19,label = 'e', lw = 3, color = 'tab:green')
    axes[0].plot(rho_pol, n_D/1e19, label = 'D', lw = 3, color = 'tab:blue')

    axes[1].plot(rho_pol, T_e, label = 'e', lw = 3, color = 'tab:green')
    axes[1].plot(rho_pol, T_D, label = 'D', lw = 3, color = 'tab:blue')

    axes[0].set_ylabel(r'n$_e$ ($10^{19}m^{-3}$)')
    axes[0].set_xlabel(r'$\rho_{pol}$')
    axes[0].grid()
    axes[0].set_xlim(right=1.025)
    #axes[0].set_yticks([0,.5,1,1.5,2,2.5])
    axes[0].legend()

    axes[1].set_ylabel(r'T$_e$ (keV)')
    axes[1].set_xlabel(r'$\rho_{pol}$')
    axes[1].grid()
    #axes[1].set_yticks([0,.5,1,1.5,2,2.5])
    axes[1].set_xlim(right=1.05)
    axes[1].legend()

    axes[0].set_ylim([0,10])
    axes[1].set_ylim([0,10])

fig.tight_layout()
plt.show()