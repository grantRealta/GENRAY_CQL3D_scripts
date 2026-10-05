"""
Plots a the ray trajectories where the colorbar corresponds to where power is absorbed.
If no power is absorbed at a given location, that part of the ray trajectory is not drawn
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.interpolate import interp2d
import matplotlib
import os, sys
from scipy.signal import find_peaks
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getGfileDict
import netCDF4
import getTargetInfo
import helperFunctions as helper
import getInputFileDictionary

gfileDict = getGfileDict.getGfileDict()
genray_in = getInputFileDictionary.getInputFileDictionary('genray_LH')
cqlinput = getInputFileDictionary.getInputFileDictionary('cql3d')

targetDir = getTargetInfo.getTargetDir()
print(targetDir)
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('legend', fontsize = 14)

#plots ray trajectories
def plotRays(toroidal = False, poloidal = True):
    #rakata
    maxDelPwrPlot = 0.95#what portion of ray power must have been damped before we stop plotting that ray

    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    rbbbs = gfileDict["rbbbs"] #R points of the LCFS
    zbbbs = gfileDict["zbbbs"] # Z points of the LCFS
    wr  = cqlrf_nc.variables["wr"][:] #major radius of the ray at each point along the trace
    wz  = cqlrf_nc.variables["wz"][:] #height of the ray at each point along the trace
    delpwr= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point
    wr *= .01; wz*=.01 #convert to m from cm

    norm = plt.Normalize(0, .03)

    figsize = (5.25,7.1)
    if machine =='WEST':
        figsize = (5.25,5)
    fig_pol,ax_pol = plt.subplots(figsize = figsize)
    plt.subplots_adjust(left=0.22,bottom = .1)
    ax_pol.set_ylabel("Z (m)")
    ax_pol.set_xlabel("R (m)")
    ax_pol.plot(xlim, ylim, 'r', lw = 2)#plot wall
    ax_pol.plot(rbbbs, zbbbs, 'k', lw = 1.5)#plot LCFS    
    
    #plot the ray using a LineCollection which allows the colormap to be applied to each ray
    for ray in range(len(wr)):
        rayPower = delpwr[ray]/delpwr[ray][0]
        powerChange = -(rayPower[1:] - rayPower[:-1])
        powerChange[powerChange < .002] = np.nan #ignore very minor changes in ray power

        mostPowerDep = helper.findNearestIndex(1 - maxDelPwrPlot, rayPower) #find the index of the last ray point we want to plot

        points = np.array([wr[ray][:mostPowerDep-1], wz[ray][:mostPowerDep-1]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # Create a continuous norm to map from data points to colors
        lc = LineCollection(segments, norm = norm,cmap=plt.cm.jet)
        # Set the values used for colormapping
        lc.set_array(powerChange[:mostPowerDep])
        lc.set_linewidth(1)
        added = ax_pol.add_collection(lc)

    #"""
    avgSPA_allLobes, onePassDelpwr, initialDelPwr = helper.getSPA(targetDir)
    avgSPA_forwardLobe = avgSPA_allLobes[0]
    
    print(f'SPA for each lobe: {avgSPA_allLobes}')
    N_para_launch = (genray_in['grill']['anmax(1)'] + genray_in['grill']['anmin(1)'])/2

    ax_pol.set_title(r'N$_{\parallel, LCFS}$ = ' + f'{N_para_launch:.2f}\n'+r'SPA$_{forward}$ = ' + f'{avgSPA_forwardLobe:.3f}')#, Shot {shotNum}')
    ax_pol.set_aspect('equal')
    if machine != 'FENIX':
        ax_pol.set_ylim(min(ylim)*1.05, max(ylim)*1.05)
        ax_pol.set_xlim(min(xlim)*.95, max(xlim)*1.05)
    else:
        ax_pol.set_ylim(-.35,.35)
        ax_pol.set_xlim(.5,.95)
    cbar = fig_pol.colorbar(added, ax = ax_pol, shrink = .9, ticks = [], pad = .01)
    cbar.set_label(r"Change in ray power (arb units)")

    helper.drawFluxSurfaces(ax_pol)
    fig_pol.tight_layout()


def main():
    plotRays(poloidal = True, toroidal = False)
    plt.show()

main()
