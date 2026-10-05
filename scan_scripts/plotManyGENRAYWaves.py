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
import netCDF4

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
def addRays(ax, targetDir, rays = [0], label = ''):
    genray_nc = netCDF4.Dataset(f'{targetDir}/genray_ec.nc','r')

    wx  = genray_nc.variables["wx"][:]/100 #major radius of the ray at each point along the trace, in m
    wz  = genray_nc.variables["wz"][:]/100 #height of the ray at each point along the trace, in m
    delpwr= genray_nc.variables["delpwr"][:] #power in the ray at each point
    maxDelPwrPlot = .99#0.7 #what portion of ray power must have been damped before we stop plotting that ray

    norm = plt.Normalize(0, 1)
    #"""

    for ray in rays:
        if len(wx[ray]) <= 1:
            continue
        mostPowerDep = findNearestIndex(1 - maxDelPwrPlot, delpwr[ray]) #find the index of the last ray point we want to plot
        ax.plot(wz[ray][:mostPowerDep], wx[ray][:mostPowerDep], lw = 2, label = label)

def setupPlot(ax, targetDir):
    gfileDict = getGfileDict.getGfileDict(targetDir = targetDir)
    genray_nc = netCDF4.Dataset(f'{targetDir}/genray_ec.nc','r')

    rlim = gfileDict["xlim"] #R points of the wall
    zlim = gfileDict["ylim"] #Z points of the wall
    
    bmodprofxz = genray_nc.variables["bmodprofxz"][:]
    eqdsk_x = genray_nc.variables["eqdsk_x"][:]
    eqdsk_z = genray_nc.variables["eqdsk_z"][:]

    freq = genray_nc.variables["freqcy"][:]
    omega = 2*np.pi*freq

    m_e = 9.109e-31
    q=1.602e-19

    ax.set_ylabel("r (m)")
    ax.set_xlabel("z (m)")

    print(f'resonant at {omega*m_e/q} T ')

    ax.contour(eqdsk_z, eqdsk_x, bmodprofxz.T, levels = [omega*m_e/q], colors = 'cyan', linewidths = 2)

    ax.plot(zlim, rlim, 'grey', lw = 2)#plot wall
    ax.plot(zlim, -rlim, 'grey', lw = 2)#plot wall

    ax.set_title(f"f$_{{EC}}$ = {freq/1e9:.0f} GHz")
    ax.set_aspect('equal')

    ax.set_xlim([.4,.9])
    ax.set_ylim([-.12,.12])

    helper.drawFluxSurfaces(ax, 
                            targetDir=targetDir,
                            rhosToPlot = [.5, 1])
    ax.axhline(0, color = 'k', lw = 2.25)

def plotManyWaves():
    #"""
    zsts = [.6,.65,.7]
    rst = 0.08
    betasts = [
        [-32,-37,-42,-48,-53],
        [-53,-57,-61,-65,-69],
        [-65,-68,-70,-72,-75],
        ]#[-65,-68,-70,-72,-75]#
    #[-32,-37,-42,-48,-53]

    fig,ax = plt.subplots(figsize = (10,6))

    for i,zst in enumerate(zsts):
        targetStem = f'/global/homes/g/grutherf/codes/IPS_scripts/WHAM_batch/WHAM_HTS/35GHz/WHAM_HTS_{zst}zst_{rst}rst'

        for betast in betasts[i][:1]:
            targetDir = f'{targetStem}_{betast}betast_1e19_40eV'
            print(targetDir)
            addRays(ax,targetDir, rays = [0], label = f'zst,rst = {zst,rst}, betast = {betast}')
    """

    zst = 0.6
    rst = 0.075
    targetStem = f'/global/homes/g/grutherf/codes/IPS_scripts/WHAM_batch/WHAM_HTS/110GHz_X2_1e19_40eV/WHAM_HTS_{zst}zst_{rst}rst'
    betasts = [28,36,43,50,57]

    fig,ax = plt.subplots(figsize = (10,6))

    for betast in betasts:
        targetDir = f'{targetStem}_{betast}betast_1e19_40eV'
        print(targetDir)
        addRays(ax,targetDir, rays = [0], label = f'{betast}')
    """
    setupPlot(ax, targetDir)

    ax.legend(loc = 'upper right')
    fig.tight_layout()
    plt.show()

plotManyWaves()
