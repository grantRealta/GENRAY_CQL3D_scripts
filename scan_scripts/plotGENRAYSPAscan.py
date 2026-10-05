import numpy as np
import matplotlib.pyplot as plt
import os, sys
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import helperFunctions as helper

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 18)
plt.rc('legend', fontsize = 14)
plt.rc('figure', titlesize = 16)
plt.rcParams['lines.linewidth'] = 2

zst = 0.67
rst = 0.15

betasts = [-46,-45,-44,-43,-42,-41,-40,-38,-37,-36]
SPAs = np.zeros(len(betasts))
fig,ax = plt.subplots()

targetStem = f'/global/homes/g/grutherf/codes/IPS_scripts/WHAM_batch/WHAM_HTS/35GHz/WHAM_HTS'

for i, betast in enumerate(betasts):
    targetDir = f'{targetStem}_{zst}zst_{rst}rst_{betast:.0f}betast_1e19_40eV'

    SPAs[i] = helper.getSPA_GEN(targetDir)

ax.plot(betasts, SPAs)

ax.set_ylabel('SGENRAY single pass absorption')
ax.set_xlabel('betast (degrees)')
ax.set_title(f'Launch at zst, rst = {zst, rst}')
ax.set_ylim([0,1.05])
fig.tight_layout()

plt.show()