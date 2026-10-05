"""
Plots the poloidal and/or toroidal ray traces
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib
from matplotlib.path import Path
import netCDF4

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
import getTargetInfo

gfileDict = getGfileDict.getGfileDict()

rmaxis = gfileDict["rmaxis"]
bcentr = gfileDict["bcentr"]
rzero = gfileDict["rzero"]

B0 = bcentr*rzero/rmaxis

print(f'B0: {B0}, bcentr: {bcentr}')
print(f'R0: {rmaxis}, rzero: {rzero}')

genray_in = getInputFileDictionary.getInputFileDictionary('genray_LH')
cqlinput = getInputFileDictionary.getInputFileDictionary('cql3d')

targetDir = getTargetInfo.getTargetDir()
print(targetDir)
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('figure', titlesize = 22)
plt.rc('legend', fontsize = 14)

#plots either the toroidal and/or poloidal ray trajectories
def plotRays(toroidal = False, poloidal = True, includeTitle = True):
    #rakata
    maxDelPwrPlot = 0.95#what portion of ray power must have been damped before we stop plotting that ray

    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    toroidalAngle = cqlrf_nc.variables["wphi"][:] 
    wr  = cqlrf_nc.variables["wr"][:] #major radius of the ray at each point along the trace
    wz  = cqlrf_nc.variables["wz"][:] #height of the ray at each point along the trace
    delpwr= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point
    wr *= .01; wz*=.01 #convert to m from cm

    radialVariable = (np.copy(cqlrf_nc.variables["spsi"]))#radial variable. I think it's rho_pol in my case. Doesn't really matter in this application
    norm = plt.Normalize(0, 1)

    fig_tor = None; fig_pol = None
    ax_tor = None; ax_pol = None
    if poloidal:
        if includeTitle:
            figsize = (5.25,6.5)
        else:
            figsize = (5,6)
        if machine =='WEST':
            figsize = (5.25,5)
        fig_pol,ax_pol = plt.subplots(figsize = figsize)
        plt.subplots_adjust(left=0.22,bottom = .1)
        ax_pol.set_ylabel("Z (m)")
        ax_pol.set_xlabel("R (m)")
        ax_pol.plot(xlim, ylim, color = 'grey', lw = 3)#plot wall
    if toroidal:
        fig_tor, ax_tor = plt.subplots(figsize = (6.5,5.5))
        ax_tor.set_ylabel("Y (m)")
        ax_tor.set_xlabel("X (m)")
        ax_tor.set_title('Toroidal ray')
        thetas = np.linspace(0,2*np.pi,100)
        ax_tor.plot(np.cos(thetas), np.sin(thetas), color = 'k', linewidth = 2)
        ax_tor.plot(np.cos(thetas)*2.355, np.sin(thetas)*2.355, color = 'k', linewidth = 2)
    
    #plot the ray using a LineCollection which allows the colormap to be applied to each ray
    for ray in range(len(wr)):
        
        # this if statement will let you control if you only want to plot one lobe
        #just change the 3 to whatever index you need
        #if not(genray_in['grill'][f'anmax(3)'] >= nparas[ray][0] >= genray_in['grill'][f'anmin(3)']):
        #    continue

        bounceIndex = helper.findBounceIndex(radialVariable[ray],bounceToFind = 1)
        mostPowerDep = helper.findNearestIndex(1 - maxDelPwrPlot, delpwr[ray]/delpwr[ray][0]) #find the index of the last ray point we want to plot
        #choose whenther you want to end plotting after the first bounce or once enough power is deposited
        endingIndex = mostPowerDep#min(bounceIndex,mostPowerDep)
        
        if poloidal:
            points = np.array([wr[ray][:endingIndex], wz[ray][:endingIndex]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            # Create a continuous norm to map from data points to colors
            lc = LineCollection(segments, norm = norm,cmap='turbo',zorder = 10)
            # Set the values used for colormapping
            lc.set_array(delpwr[ray][:endingIndex]/delpwr[ray][0])
            lc.set_linewidth(2)
            ax_pol.add_collection(lc)

        if toroidal:
            rayX = (np.cos(toroidalAngle[ray])*wr[ray])[:endingIndex]
            rayY = (np.sin(toroidalAngle[ray])*wr[ray])[:endingIndex]
            points = np.array([rayX, rayY]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            # Create a continuous norm to map from data points to colors
            lc = LineCollection(segments, norm = norm,cmap=plt.cm.turbo)
            # Set the values used for colormapping
            lc.set_array(delpwr[ray][:endingIndex]/delpwr[ray][0])
            lc.set_linewidth(1)
            ax_tor.add_collection(lc)

    avgSPA_allLobes, onePassDelpwr, initialDelPwr = helper.getSPA(targetDir)
    
    print(f'SPA for each lobe: {avgSPA_allLobes}')
    totalSPA = 1- np.sum(onePassDelpwr)/np.sum(initialDelPwr)
    print(f'avg SPA across all lobes: {totalSPA}')

    N_para_launch = (genray_in['grill']['anmax(1)'] + genray_in['grill']['anmin(1)'])/2

    cmap = matplotlib.cm.ScalarMappable(norm = matplotlib.colors.Normalize(0,1),
            cmap = plt.get_cmap('turbo'))
    cmap.set_array([])
    cticks = np.linspace(0,1,5)

    if poloidal:
        if includeTitle:
            fig_pol.suptitle(r'N$_{\parallel, LCFS}$ = ' + f'{N_para_launch:.2f}\n'+r'SPA$_{forward}$ = ' + f'{avgSPA_allLobes[0]:.3f}')#, Shot {shotNum}')
        ax_pol.set_aspect('equal')
        if machine != 'FENIX':
            ax_pol.set_ylim(min(ylim)-.05, max(ylim)+.05)
            ax_pol.set_xlim(min(xlim)-.05, max(xlim)+.05)
        else:
            ax_pol.set_ylim(-.35,.35)
            ax_pol.set_xlim(.5,.95)
        cbar = fig_pol.colorbar(cmap, ax = ax_pol, shrink = .8, ticks = cticks, pad = .01)
        cbar.set_label(r"Fractional power remaining in ray")

        lim = np.array([xlim, ylim]).T
        limiterPath  = Path(lim)
        
        helper.drawFluxSurfaces(ax_pol, colors = 'k', limPath = limiterPath)
        fig_pol.tight_layout()

    if toroidal:
        ax_tor.set_title(r'N$_{\parallel, launch}$ = ' + f'{N_para_launch:.1f}')
        ax_tor.set_aspect('equal')
        thetas = np.linspace(0,2*np.pi,100)
        ax_tor.plot(np.cos(thetas), np.sin(thetas), color = 'k', linewidth = 2)
        ax_tor.plot(np.cos(thetas)*2.355, np.sin(thetas)*2.355, color = 'k', linewidth = 2)
        cbar = fig_tor.colorbar(cmap, ax = ax_tor, shrink = .75, ticks = cticks, pad = .01)
        cbar.set_label(r"Fractional power remaining in ray")
        fig_tor.tight_layout()

#plots where power is deposited in the first pass
#it's a bit rough, but it gives a good sense
#this is also probably deprecated
def plotFirstPassAbsorption():
    radialBinEdges = np.linspace(0,1,51)
    radialBinCenters = (radialBinEdges[1:]+radialBinEdges[:-1])/2
    powerDep = np.zeros(len(radialBinCenters))
    delpwr= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point
    radialVariable = (np.copy(cqlrf_nc.variables["spsi"]))
    fig, ax = plt.subplots()
    nparas = cqlrf_nc.variables['wnpar'][:]
    for ray in range(len(nparas)):
        if genray_in['grill']['anmax(1)'] >= nparas[ray][0] >= genray_in['grill']['anmin(1)']:
            mostPowerDep = helper.findNearestIndex(1 - 0.99, delpwr[ray]/delpwr[ray][0]) #find the index of the last ray point we want to plot
            firstBounceIndex = helper.findBounceIndex(radialVariable[ray][:mostPowerDep],bounceToFind = 1)

            radialVariableCenters = (radialVariable[ray][:firstBounceIndex][1:] + radialVariable[ray][:firstBounceIndex][:-1])/2

            indices = np.digitize(radialVariableCenters, radialBinEdges, right = False)
            indices[indices>49] = 49
            powerDep[indices-1] += np.diff( delpwr[ray][:firstBounceIndex])
    
    ax.plot(radialBinCenters, powerDep)
    ax.set_xlim([0,1])


def main():
    plotRays(poloidal = True, toroidal = True, includeTitle = False)
    #plotFirstPassAbsorption()
    plt.show()

main()
