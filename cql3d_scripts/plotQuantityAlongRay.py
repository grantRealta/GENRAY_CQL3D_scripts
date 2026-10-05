"""
Useful script for plotting a given quantity along a ray
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import netCDF4
import helperFunctions as helper
import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')
genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')

import getGfileDict
gfileDict = getGfileDict.getGfileDict()

plt.rc('xtick', labelsize = 17)
plt.rc('ytick', labelsize = 17)
plt.rc('axes', labelsize = 19)
plt.rc('axes', titlesize = 18)
plt.rc('legend', fontsize = 16)

rya = cql_nc.variables["rya"][:]
q_prof = cql_nc.variables["qsafety"][:]
qFunc = interp1d(rya, q_prof,fill_value = 'extrapolate', bounds_error = False)

Nparas = np.copy(genray_nc.variables["wnpar"])
Nphis = np.copy(genray_nc.variables["wn_phi"])
poloidalDistance = np.copy(genray_nc.variables["ws"]) 
radialCoord = np.copy(genray_nc.variables["spsi"]) 
sbtot = np.copy(genray_nc.variables["sbtot"]) /1e4 #convert to T
sbr = np.copy(genray_nc.variables["sb_r"]) /1e4 #convert to T
sbz = np.copy(genray_nc.variables["sb_z"]) /1e4 #convert to T
sbpol = -np.sqrt(sbr**2 + sbz**2)
sene =  np.copy(genray_nc.variables["sene"]) *1e6#convert to m^-3
delpwr= np.copy(cqlrf_nc.variables["delpwr"]) #power in the ray at each point
thetas = np.copy(genray_nc.variables["w_theta_pol"])*np.pi/180
Rs = np.copy(genray_nc.variables["wr"])/1e2#cm to m
Zs = np.copy(genray_nc.variables["wz"])/1e2
eps_0=8.85*10**-12
q = 1.602e-19
m_e = 9.109e-31
m_D = 3.343e-27/2
w_pis = np.sqrt(sene*q**2/(eps_0*m_D))
w_pes = np.sqrt(sene*q**2/(eps_0*m_e))
W_is = q*sbtot/m_D
W_es = q*sbtot/m_e
c = 3e8
#w_LHs = 1/np.sqrt(1/w_pis**2 + 1/(W_is * W_es))
w = 2*np.pi*4.6e9

Rmaxis = gfileDict['rmaxis']
Zmaxis = gfileDict['zmaxis']
Zs_rel = Zs - Zmaxis
Rs_rel = Rs - Rmaxis

rs = np.sqrt(Zs_rel**2 + Rs_rel**2)

fig,ax = plt.subplots()

minRatioToPlot = 0.8
rayOfInterest = [18]
for i in rayOfInterest:#range(len(Nparas)):
    Npara = Nparas[i]
    kpara = Npara*(w/c)
    Nphi = Nphis[i]
    kphi = Nphi*(w/c)
    thetaComp = (Npara-Nphi)*(w/c)
    r = rs[i]
    freqTerm = (w_pes**2/(w_pes**2 + W_es**2))[i]

    mdot = r*(w*np.sin(thetas[i])/Rs[i])*(freqTerm + (1/kpara)*(thetaComp - kphi))
    ray = i
    delpwrRatios = delpwr[ray]/np.max(delpwr[ray])
    endIndex = helper.findNearestIndex(minRatioToPlot, delpwrRatios) 
    #ax.plot(poloidalDistance[ray][:endIndex], sbtot[ray][:endIndex])
    ax.plot(radialCoord[ray][:endIndex], mdot[:endIndex], lw = 2)
print(sbpol[18])
ax.set_title(f'Shot {shotNum}')
ax.set_ylabel(r'$B_\theta$')
ax.set_ylabel(r'$\dot{m}$')
ax.set_xlabel(r'$\rho_{pol}$')
#ax.set_xlabel('Poloidal distance along ray (cm)')
#ax.set_ylim([1,6])
ax.set_xlim([0,1])
#ax.set_ylim([-1.5e8,4.5e8])
fig.tight_layout()
plt.show()