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

#returns the index of the array whose element is closest to value
def findNearestIndex(value, array):
    idx = (np.abs(array - value)).argmin()

    return idx

#adds the ray traces to ax
def plotStix():

    toPlot = 'Nperp2_plus'

    plotRoot = True

    genrayPath = f'{targetDir}/genray.in'
    nml = f90nml.read(genrayPath)

    betast = nml['eccone']['betast']


    rlim = gfileDict["xlim"] #R points of the wall
    zlim = gfileDict["ylim"] #Z points of the wall
    
    bmodprofxz = genray_nc.variables["bmodprofxz"][:]
    densprofxz = genray_nc.variables['densprofxz'][:] * 1e6 #m-3

    print(f'max of bmodprofxz and densprofxz: {np.max(bmodprofxz), np.max(densprofxz)}')

    eqdsk_x = genray_nc.variables["eqdsk_x"][:]
    eqdsk_z = genray_nc.variables["eqdsk_z"][:]
    Z,X = np.meshgrid(eqdsk_x,eqdsk_z)

    freq = genray_nc.variables["freqcy"][:]
    omega = 2*np.pi*freq

    m_e = 9.109e-31
    m_D = 3.343e-27
    q_e=-1.602e-19
    q_D=1.602e-19
    eps_0 = 8.854e-12

    W_ce = q_e*bmodprofxz/m_e
    W_cD = q_D*bmodprofxz/m_D
    w_pe2 = q_e**2*densprofxz/(m_e*eps_0)
    w_pD2 = q_D**2*densprofxz/(m_D*eps_0)

    S = 1 - (w_pe2/(omega**2 - W_ce**2)) - (w_pD2/(omega**2 - W_cD**2))
    D = (W_ce/omega)*(w_pe2/(omega**2-W_ce**2)) + (W_cD/omega)*(w_pD2/(omega**2-W_cD**2))
    P = 1 - (w_pe2/omega**2) - (w_pD2/omega**2)

    R = S + D
    L = S - D

    Npar = 1.25

    C_0 = P*((Npar**2-S)**2-D**2)
    C_2 = (Npar**2 - S)*(P+S) + D**2
    C_4 = S

    Nperp2_minus = (-C_2 - np.sqrt(C_2**2-4*C_4*C_0))/(2*C_4)
    Nperp2_plus = (-C_2 + np.sqrt(C_2**2-4*C_4*C_0))/(2*C_4)

    Nperp2_plus[Nperp2_plus < 0] = np.nan
    Nperp2_minus[Nperp2_minus < 0] = np.nan

    norm = plt.Normalize(0, 1)
    #"""
    fig,ax = plt.subplots(figsize = (10,6))
    plt.subplots_adjust(left=0.22,bottom = .1)
    ax.set_ylabel("r (m)")
    ax.set_xlabel("z (m)")

    print(f'resonant at {np.abs(omega*m_e/q_e)} T ')

    ax.contour(eqdsk_z, eqdsk_x, bmodprofxz.T, levels = [np.abs(omega*m_e/q_e)], colors = 'r', linewidths = 2)

    if toPlot == 'Nperp2_plus':
        pcm = ax.pcolormesh(eqdsk_z, eqdsk_x, Nperp2_plus.T, cmap = 'viridis', 
                            vmin = 0, vmax = 3)
        cb = fig.colorbar(pcm, ax = ax, shrink = 0.6)

    if toPlot == 'Nperp2_minus':
        pcm = ax.pcolormesh(eqdsk_z, eqdsk_x, Nperp2_minus.T, cmap = 'viridis', 
                            vmin = 0, vmax = 3)
        cb = fig.colorbar(pcm, ax = ax, shrink = 0.6)

    if toPlot == 'L':
        pcm = ax.pcolormesh(eqdsk_z, eqdsk_x, L.T, cmap = 'bwr', 
                            vmin = -np.nanmax(L), vmax = np.nanmax(L))
        cb = fig.colorbar(pcm, ax = ax, shrink = 0.6)

    if toPlot == 'R':
        pcm = ax.pcolormesh(eqdsk_z, eqdsk_x, R.T, cmap = 'bwr', 
                            vmin = -10, vmax = 10)
        cb = fig.colorbar(pcm, ax = ax, shrink = 0.6)

    if toPlot == 'D':
        pcm = ax.pcolormesh(eqdsk_z, eqdsk_x, D.T, cmap = 'bwr', 
                            vmin = -np.nanmax(D), vmax = np.nanmax(D))
        cb = fig.colorbar(pcm, ax = ax, shrink = 0.6)

    if toPlot == 'S':
        print(f'min, max: {np.nanmin(S), np.nanmax(S)}')
        pcm = ax.pcolormesh(eqdsk_z, eqdsk_x, S.T, cmap = 'bwr', 
                            vmin = -10, vmax = 10)
        cb = fig.colorbar(pcm, ax = ax, shrink = 0.6)

    if toPlot == 'P':
        P[P < 0] = np.nan
        pcm = ax.pcolormesh(eqdsk_z, eqdsk_x, P.T, cmap = 'viridis', 
                            vmin = 0, vmax = 1)
        cb = fig.colorbar(pcm, ax = ax, shrink = 0.6)

    if toPlot == 'RL/S':
        pcm = ax.pcolormesh(eqdsk_z, eqdsk_x, (R*L/S).T, cmap = 'bwr', 
                            vmin = -5, vmax = 5)
        cb = fig.colorbar(pcm, ax = ax, shrink = 0.6)

    ax.plot(zlim, rlim, 'grey', lw = 2)#plot wall
    ax.plot(zlim, -rlim, 'grey', lw = 2)#plot wall

    ax.set_title(f"Plotting {toPlot}, f$_{{EC}}$ = {freq/1e9:.0f} GHz, betast = {betast:.0f} deg")
    ax.set_aspect('equal')

    if machine == 'WHAM':
        pass
        ax.set_xlim([.4,.9])
        ax.set_ylim([-.12,.12])

    helper.drawFluxSurfaces(ax, 
                            targetDir=getTargetInfo.getTopmostShotDir(),
                            rhosToPlot = [.5, 1])
    ax.axhline(0, color = 'k', lw = 2.25)

    
    wr  = genray_nc.variables["wr"][:]/100 #major radius of the ray at each point along the trace, in m
    wx  = genray_nc.variables["wx"][:]/100 #major radius of the ray at each point along the trace, in m
    wz  = genray_nc.variables["wz"][:]/100 #height of the ray at each point along the trace, in m
    delpwr= genray_nc.variables["delpwr"][:] #power in the ray at each point
    maxDelPwrPlot = 1#0.7 #what portion of ray power must have been damped before we stop plotting that ray

    norm = plt.Normalize(0, 1)

    #"""
    for ray in range(len(wr)):
        if len(wr[ray]) <= 1:
            continue

        delpwr[ray,:] = delpwr[ray,:]/delpwr[ray,0] #normalize the ray power to that ray's starting power
        mostPowerDep = findNearestIndex(1 - maxDelPwrPlot, delpwr[ray]) #find the index of the last ray point we want to plot

        points = np.array([wz[ray][:mostPowerDep], wx[ray][:mostPowerDep]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Create a continuous norm to map from data points to colors
        lc = LineCollection(segments, norm = norm,cmap=plt.cm.turbo)
        # Set the values used for colormapping
        lc.set_array(delpwr[ray][:mostPowerDep])
        lc.set_linewidth(3)
        ax.add_collection(lc)


    fig.tight_layout()
    plt.show()

plotStix()
