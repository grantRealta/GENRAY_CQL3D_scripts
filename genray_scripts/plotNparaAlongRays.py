
"""
Plots either the N|| evolution of all rays (with the color bar corresponding to the power remaining in the rays)
of plots the N|| evolution of a single ray, as well as the poloidal and toroidal contribution to N|| and Nacc
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.optimize import fsolve

print('past imports')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 14)

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import helperFunctions as helper
import netCDF4

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
genray_nc = netCDF4.Dataset(f'{targetDir}/genray_ec.nc','r')

def plotNEvolution():

    m_e = 9.109e-31
    m_D = 3.343e-27
    q=1.602e-19
    eps_0 = 8.854e-12

    Nparas = np.copy(genray_nc.variables["wnpar"]) #n_|| of the ray at each point along the ray trace
    Nperps = np.copy(genray_nc.variables["wnper"]) #n_perp of the ray at each point along the ray trace
    sbtot = np.copy(genray_nc.variables["sbtot"])/1e4 #B along ray
    sene = np.copy(genray_nc.variables["sene"])*1e6 #ne along ray
    delpwr = np.copy(genray_nc.variables["delpwr"]) #ne along ray

    swpe = np.sqrt(sene*q**2/(eps_0*m_e))

    ws= np.copy(genray_nc.variables["ws"])/1e2 #distance along the ray

    fig,ax = plt.subplots()#,figsize = (7,7), sharex = True)

    for j in range(len(ws)):
        if j > 0:
            pass
        normedDelpwr = delpwr[j]
        normedDelpwr = normedDelpwr/np.nanmax(normedDelpwr)
        #ax.plot(ws[j][ws[j] > 0], Nparas[j][ws[j] > 0])
        ax.plot(1-normedDelpwr[normedDelpwr > 0.01], Nparas[j][normedDelpwr > 0.01])
        #ax.plot(ws[j][ws[j] > 0], swpe[j][ws[j] > 0], lw = 2)

    #ax.axhline(2*np.pi*28e9, color = 'k', linestyle = 'dashed', lw = 2)

    ax.set_ylabel('N||')
    #ax.set_xlabel('distance along ray (m)')
    ax.set_xlabel('Fraction of power damped')
    #ax.set_ylim([-1.3,.5])
    fig.tight_layout()

    plt.show()

plotNEvolution()
