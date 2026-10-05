###
# 2D plot of B
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
from matplotlib.lines import Line2D


plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('legend', fontsize = 14)
plt.rc('figure', titlesize = 16)

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

xlim = gfileDict["xlim"] #R points of the wall
ylim = gfileDict["ylim"] #Z points of the wall
rbbbs = gfileDict["rbbbs"] #R points of the LCFS
zbbbs = gfileDict["zbbbs"] # Z points of the LCFS

B_rGrid[np.isnan(B_rGrid)] = 0
B_TGrid[np.isnan(B_TGrid)] = 0
B_zGrid[:,0] = B_zGrid[:,1]

B_totstrength = np.sqrt(np.square(B_zGrid) + np.square(B_rGrid) + np.square(B_TGrid))

limFunc = interp1d(ylim, xlim)
Zs, Rs = np.meshgrid(zgrid, rgrid)
Rs_interp = limFunc(Zs)
B_totstrength.T[Rs_interp < Rs] = np.nan

rgrid = np.concatenate((-1*np.flip(rgrid[1:]), rgrid))

# Mirror B in the z direction
B_totstrength = np.concatenate(
    (np.flip(B_totstrength[:, 1:], axis=1), B_totstrength),
    axis=1
)

print(f'min, max of B: {np.nanmin(B_totstrength), np.nanmax(B_totstrength)}')

on_axisB = B_totstrength[int(len(zgrid)/2), int(len(rgrid)/2)]

print(f'on_axisB: {on_axisB}')


fig,ax = plt.subplots(figsize = (10,4))

#p1 = ax.pcolormesh(zgrid, rgrid, B_totstrength.T,shading = 'nearest',cmap='plasma')
p1 = ax.pcolormesh(zgrid, rgrid, B_totstrength.T/on_axisB,shading = 'nearest',cmap='plasma',vmax = 20)

plotRes = True

if plotRes:
    q = 1.602e-19
    m_H = 1.67e-27
    m_e = 9.109e-31

    ICRH = False
    if ICRH:
        fs = np.array([35,50,65])
        factor = 1e6
        B_res = fs*factor*np.pi*2*m_H/q
        colors = ['r', 'cyan', 'lime', 'white', 'black', 'magenta']
        suffix = 'MHz'

    else:
        fs = np.array([170])#np.array([35,110])
        factor = 1e9
        colors = ['r', 'cyan', 'lime', 'white', 'black', 'magenta']
        B_res = fs*factor*np.pi*2*m_e/q
        suffix = 'GHz'

ax.contour(zgrid, rgrid, B_totstrength.T, levels = B_res, lw = 2, colors = colors)


handles = [
    Line2D([0], [0], color=color, lw=2.25, label=f"{f} {suffix}")
    for color, f in zip(colors, fs)
]

ax.set_ylim([-.2,.2])
ax.set_xlim([.5,.8])

ax.set_aspect('equal')
ax.fill_between(ylim, xlim, np.ones(len(ylim))*np.max(rgrid),color='w')#plot wall
ax.fill_between(ylim, -1*xlim, -1*np.ones(len(ylim))*np.max(rgrid),color='w')#plot wall
ax.plot(zbbbs, rbbbs, 'k', lw = 1.5)#plot LCFS
ax.plot(zbbbs, -rbbbs, 'k', lw = 1.5)#plot LCFS
cbar = fig.colorbar(p1, ax = ax, shrink = .7, pad = .01)
#cbar.set_label(r"B$_{total}$ (T)")
cbar.set_label(r"B$_{min}$/B$_{min}$ (T)")
ax.set_ylabel('r (m)', labelpad = -10)
ax.set_xlabel('z (m)')
if machine == 'WHAM':
    ax.set_xlim([-1.1,1.1])
else:
    ax.set_xlim([np.min(zgrid),np.max(zgrid)])

ax.set_ylim([-np.max(rgrid),np.max(rgrid)])


ax.legend(handles = handles, loc = 'upper right', ncol = 2)
fig.tight_layout()
plt.show()
