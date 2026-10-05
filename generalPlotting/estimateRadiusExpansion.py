###
# Estimates the change in radius at finite beta
###

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import getGfileDict
import getTargetInfo


plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 18)
plt.rc('legend', fontsize = 14)
plt.rc('figure', titlesize = 16)
plt.rcParams['lines.linewidth'] = 2

targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()
print(targetDir)

gfileDict = getGfileDict.getGfileDict()

#get the variables we'll need for plotting
rgrid = gfileDict["rgrid"]
zgrid = gfileDict["zgrid"]
B_zGrid = gfileDict["bzrz"]
B_TGrid = gfileDict["btrz"]
B_rGrid = gfileDict["brrz"]

B_rGrid[np.isnan(B_rGrid)] = 0
B_TGrid[np.isnan(B_TGrid)] = 0
B_zGrid[:,0] = B_zGrid[:,1]

rbbbs = gfileDict["rbbbs"] #R points of the LCFS
zbbbs = gfileDict["zbbbs"] # Z points of the LCFS
limFunc = interp1d(zbbbs, rbbbs)

#this is ordered (z, r)
B_totstrength = np.sqrt(np.square(B_zGrid) + np.square(B_rGrid) + np.square(B_TGrid))


limRs = limFunc(zgrid)
onAxisB = B_totstrength[:, 0]

Bfunc = interp1d(zgrid, onAxisB)

fig,ax = plt.subplots()
ax.plot(rgrid, B_totstrength[int(len(zgrid)/2), :])
ax.set_ylabel('|B| (T)')
ax.set_xlabel('R (m)')
fig.tight_layout()
plt.show()

B_midplane = onAxisB[int(len(onAxisB)/2)]

if machine == 'Hammir':
    z_L = 0
    r_L = limFunc(z_L)#.22527
    B_L = Bfunc(z_L)

    z_th = 5
    r_th = limFunc(z_th)
    B_Th = Bfunc(z_th)

R_m = B_Th/onAxisB
#loss cone angles
thetas = np.arcsin(1/np.sqrt(R_m))

beta_prof_normed = ((np.pi - thetas)/np.pi) * (1/onAxisB**2)
beta_prof_normed[np.isnan(beta_prof_normed)] = 0

fig,ax = plt.subplots()
betas = [0,.2,.4]

r_zeroBeta = r_L*np.sqrt(B_L/onAxisB)


for beta in betas:
    #newr = r_th*np.sqrt(R_m)*(1-beta)**(-.25)
    beta_prof = beta_prof_normed*beta/beta_prof_normed[int(len(beta_prof_normed)/2)]
    print(f'max of beta prof: {np.nanmax(beta_prof)}')
    beta_L = interp1d(zgrid, beta_prof)(z_L)
    print(f'beta_L: {beta_L}')



    newr = r_L*np.sqrt(B_L/onAxisB)*((1-beta_L)/(1-beta_prof))**.25
    #print(f'max of new r: {np.nanmax(newr)}')
    #newr = limRs/((1-beta_prof)**(.25))
    ax.plot(zgrid, newr-r_zeroBeta, label = r'$\beta = $'+f'{beta}')
#ax.plot(zgrid, limRs, label = 'basecase', color = 'k')


ax.legend(loc = 'best')
#ax.set_ylabel(r'separatrix radius (m)')
ax.set_ylabel(r'$\Delta r = r(\beta) - r(\beta = 0)$ (m)')
ax.set_xlabel('z (m)')
ax.set_xlim([-5,5])

ax.set_title(f'limiter at z = {z_L:.2f} m, r = {r_L:.3f} m')

fig.tight_layout()
plt.show()