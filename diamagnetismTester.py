import matplotlib.pyplot as plt
import numpy as np

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('figure', titlesize = 22)
plt.rc('legend', fontsize = 14)

from pleiades.configurations import WHAM
from pleiades import Gaussian, Flat, KesnerKotelnikov
print('past imports')
#initialize wham configuration
mu0 = 1.256637e-6
beta0 = 0.4
wham = WHAM()

wham.vacuum_fields()
wham.working_mesh = wham.PlasmaMesh

B_vac = wham.B
nmid = B_vac.shape[0]//2
psi_vac = wham.psi_n
B_0 = B_vac[nmid,0]
p_0 = B_0**2/(2*mu0)*beta0
psi_vac = np.copy(wham.psi)
psilim_vac = wham.psilim

wham.plasma.radial_profile = Flat()
wham.plasma.axial_profile = KesnerKotelnikov()
wham.plasma.p_max = p_0
wham.plasma.p_kin = 0.0
wham.solve_equil()

psi_pla = np.copy(wham.psi)
lvls1 = np.linspace(0.0,np.sqrt(psilim_vac),7)
lvls2 = np.linspace(0.0,np.sqrt(wham.psilim),7)

fig, ax = plt.subplots()

Z = wham.Z
R = wham.R

ax.contour(wham.Z,wham.R,np.sqrt(psi_vac),levels=lvls1, colors='b')
ax.contour(wham.Z,wham.R,np.sqrt(psi_pla),levels=lvls2, colors='r')

ax.contour(wham.Z,-1*wham.R,np.sqrt(psi_vac),levels=lvls1, colors='b')
ax.contour(wham.Z,-1*wham.R,np.sqrt(psi_pla),levels=lvls2, colors='r')

ax.plot([0],[0], label = r'$\beta = 0$ %', color = 'b')
ax.plot([0],[0], label = r'$\beta = 40$ %', color = 'r')

ax.set_ylim([-.15,.15])
ax.legend()
ax.set_ylabel('r (m)')
ax.set_xlabel('z (m)')
fig.tight_layout()

plt.show()