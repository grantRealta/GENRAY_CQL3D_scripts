"""
This is effectively equivalent to plotDistributionFunction_2D.py, but it plots slices are desired radii
"""

import numpy as np
import matplotlib.pyplot as plt
import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import netCDF4
import getTargetInfo
import helperFunctions as helper
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

shotNumber = shotNum.split('.')[0]
shotTime = shotNum.split('.')[1][1:]

print(targetDir)

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')

plt.rc('xtick', labelsize = 12)
plt.rc('ytick', labelsize = 12)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 13)

def plotEnergeticHeatMap():
    rya = cql_nc.variables["rya"][:]#radial points at which f was solved

    #distribution function (lrz, jx, iy)
    #vnorm^3*s^3/cm^6
    f = cql_nc.variables["f"][:]

    #pitch angles mesh at which f is defined
    pitchAngleMesh = np.ma.getdata(cql_nc.variables["y"][:])
    #normalized speed mesh of f
    normalizedVel = cql_nc.variables["x"][:]

    #energy  = restmkev*(gamma-1)
    #energies corresponding to velocities jx
    enerkev = cql_nc.variables["enerkev"][:] 

    #minimum energy of particles to consider for plotting
    minEnergy = 50
    #index of that minimum energy in enerkev
    #this index is also the index for the corresponding velocity
    minEnergyIndex = 0#np.where(enerkev < minEnergy)[0][-1]
    #distribution function for energetic particles
    energeticF = f[:,:,:]

    energeticF_integOverPitch = np.zeros((len(rya), len(enerkev[minEnergyIndex:])))
    for rhoIndex in range(len(rya)):
        #this is the angular part of the spherical jacobian
        integOverPitch = 2*np.pi*np.trapz(energeticF[rhoIndex,:]*np.sin(pitchAngleMesh[rhoIndex]), pitchAngleMesh[rhoIndex], axis = 1)
        energeticF_integOverPitch[rhoIndex,:] = integOverPitch
    
    fig,ax = plt.subplots(figsize=(8,5.5))

    rhosToPlot = np.array([.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
    for rho_pol in rhosToPlot:
        rho_index = helper.findNearestIndex(rho_pol, rya)

        ax.plot(enerkev, energeticF_integOverPitch[rho_index,:], lw = 2, label = r'$\rho_{pol} = $' + f'{rho_pol}')

    ax.set_title(f'Shot {shotNumber} at {shotTime} ms ')

    ax.legend(loc='best',ncol=2)
    ax.set_ylabel(r'Distribution function (vnorm$^3\cdot$s$^3$/cm$^6$)')
    ax.set_xlabel('Energy (keV)')
    ax.set_yscale('log')
    ax.set_xlim([0,250])
    ax.set_ylim([1e10,1e14])
    fig.tight_layout()
    plt.show()

plotEnergeticHeatMap()