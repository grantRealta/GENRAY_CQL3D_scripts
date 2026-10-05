###
# Plot the ray trajectories and damping as predicted by GENRAY
# since it's a linear code, the damping is very wrong for LH, but the ray trajectory can be useful
# may also take a long time to plot since the rays are apt to be much longer than they would be when proper damping is accounted for
###

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os, sys

#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getGfileDict
import helperFunctions as helper
import getInputFileDictionary
import netCDF4
import f90nml
import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
machine = getTargetInfo.getMachine()
genray_in = getInputFileDictionary.getInputFileDictionary('genray')
gfileDict = getGfileDict.getGfileDict(targetDir = getTargetInfo.getTopmostShotDir())
genray_nc = netCDF4.Dataset(f'{targetDir}/genray_ec.nc','r')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 20)
plt.rc('figure', titlesize = 20)
plt.rc('legend', fontsize = 14)
plt.rcParams['lines.linewidth'] = 2

#adds the ray traces to ax
def plotStix():
    sene =genray_nc.variables["sene"][:]*1e6 #density along ray
    sbtot = genray_nc.variables["sbtot"][:]/1e4 #magnetic field along ray

    freq = genray_nc.variables["freqcy"][:]

    Npar = genray_nc.variables["wnpar"][:] #N||
    #Npar = np.zeros(Npar.shape)+1.25

    omega = 2*np.pi*freq

    m_e = 9.109e-31
    m_D = 3.343e-27
    q_e=-1.602e-19
    q_D=1.602e-19
    eps_0 = 8.854e-12

    W_ce = q_e*sbtot/m_e
    W_cD = q_D*sbtot/m_D
    w_pe2 = q_e**2*sene/(m_e*eps_0)
    w_pD2 = q_D**2*sene/(m_D*eps_0)

    S = 1 - (w_pe2/(omega**2 - W_ce**2)) - (w_pD2/(omega**2 - W_cD**2))
    D = (W_ce/omega)*(w_pe2/(omega**2-W_ce**2)) + (W_cD/omega)*(w_pD2/(omega**2-W_cD**2))
    P = 1 - (w_pe2/omega**2) - (w_pD2/omega**2)

    """
    cweps11 = genray_nc.variables["cweps11"][:]
    S_gen = cweps11[0] + 1j*cweps11[1]
    cweps12 = genray_nc.variables["cweps12"][:]
    D_gen = cweps12[0] + 1j*cweps12[1]
    cweps33 = genray_nc.variables["cweps33"][:]
    P_gen = cweps33[0] + 1j*cweps33[1]

    S_gen = np.real(S_gen)
    D_gen = np.imag(D_gen)
    P_gen = np.real(P_gen)
    """

    C_0 = P*((Npar**2-S)**2-D**2)
    C_2 = (Npar**2 - S)*(P+S) + D**2
    C_4 = S

    """
    C_0_gen = P_gen*((Npar**2-S_gen)**2-D_gen**2)
    C_2_gen = (Npar**2 - S_gen)*(P_gen+S_gen) + D_gen**2
    C_4_gen = S_gen
    """

    Nperp2_fast = (-C_2 - np.sqrt(C_2**2-4*C_4*C_0))/(2*C_4)
    Nperp2_slow = (-C_2 + np.sqrt(C_2**2-4*C_4*C_0))/(2*C_4)

    #Nperp2_fast_gen = (-C_2_gen - np.sqrt(C_2_gen**2-4*C_4_gen*C_0_gen))/(2*C_4_gen)
    #Nperp2_slow_gen = (-C_2_gen + np.sqrt(C_2_gen**2-4*C_4_gen*C_0_gen))/(2*C_4_gen)

    Nperp2_ray =  genray_nc.variables["wnper"][:]**2

    cpdr_residual = (C_4 * Nperp2_ray**2 + C_2 * Nperp2_ray + C_0)
    #cpdr_residual_gen = (C_4_gen * Nperp2_ray**2 + C_2_gen * Nperp2_ray + C_0_gen)
    
    ws  = genray_nc.variables["ws"][:]/100
    wr  = genray_nc.variables["wr"][:]/100

    fig, ax = plt.subplots()
    #ax.set_ylim([-2,2])
    ax.axhline(0, color = 'k', linestyle = 'dashdot')
    for ray in range(len(wr)):

        if len(wr[ray]) <= 1:
            continue

        if ray > 0:
            continue

        """
        ax.plot(ws[ray][ws[ray] > 0], C_0[ray][ws[ray] > 0] / np.max(np.abs(C_0[ray][ws[ray] > 0])))
        ax.plot(ws[ray][ws[ray] > 0], C_2[ray][ws[ray] > 0]/ np.max(np.abs(C_2[ray][ws[ray] > 0])))
        ax.plot(ws[ray][ws[ray] > 0], C_4[ray][ws[ray] > 0]/ np.max(np.abs(C_4[ray][ws[ray] > 0])))
        ax.plot(ws[ray][ws[ray] > 0], Nperp2[ray][ws[ray] > 0] / np.max(np.abs(Nperp2[ray][ws[ray] > 0])))
        """

        #ax.plot(ws[ray][ws[ray] > 0], P_gen[ray][ws[ray] > 0], color = 'purple', linestyle = 'dotted', zorder = 1000, label = 'GENRAY ')


        ax.plot(ws[ray][ws[ray] > 0], Nperp2_ray[ray][ws[ray] > 0], color = 'r', label = 'GENRAY ')


        ax.plot(ws[ray][ws[ray] > 0], Nperp2_fast[ray][ws[ray] > 0], label = 'minus root', color = 'g')
       # ax.plot(ws[ray][ws[ray] > 0], Nperp2_fast_gen[ray][ws[ray] > 0], color = 'forestgreen', linestyle = 'dashed')

        ax.plot(ws[ray][ws[ray] > 0], Nperp2_slow[ray][ws[ray] > 0],label = 'plus root', color = 'b', zorder = 100, linestyle = 'dashed')
        #ax.plot(ws[ray][ws[ray] > 0], 1+Nperp2_slow_gen[ray][ws[ray] > 0], color = 'b', linestyle = 'dashed')


        #ax.plot(ws[ray][ws[ray] > 0], cpdr_residual[ray][ws[ray] > 0], color = 'k', linestyle = 'dotted')
        #ax.plot(ws[ray][ws[ray] > 0], cpdr_residual_gen[ray][ws[ray] > 0], color = 'grey', linestyle = 'dashed')

        #ax.plot(ws[ray][ws[ray] > 0], sene[ray][ws[ray] > 0]/1e19)
        
    ax.set_ylabel('Nperp^2')
    ax.set_xlabel('distance along ray (m)')
    ax.legend(loc = 'best')


    fig.tight_layout()
    plt.show()

plotStix()
