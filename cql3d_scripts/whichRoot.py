"""
Plots whether the ray(s) are on the fast or slow branch
Options are to plot the N|| evolution and/or a poloidal ray trace
ONLY WORKS IF ID = 2
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import patches
import warnings
warnings.filterwarnings("ignore")

plt.rc('xtick', labelsize = 17)
plt.rc('ytick', labelsize = 17)
plt.rc('axes', labelsize = 19)
plt.rc('axes', titlesize = 18)
plt.rc('legend', fontsize = 16)

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getGfileDict
gfileDict = getGfileDict.getGfileDict()

import getInputFileDictionary
genInput = getInputFileDictionary.getInputFileDictionary('genray_')

import helperFunctions as helper
import netCDF4

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')
genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')

plt.rc('xtick', labelsize = 17)
plt.rc('ytick', labelsize = 17)
plt.rc('axes', labelsize = 19)
plt.rc('axes', titlesize = 18)
plt.rc('legend', fontsize = 16)

m_e = 9.109e-31
m_D = 3.343e-27
eps_0 = 8.854e-12
e = 1.6e-19
c = 3e8
w = 2*np.pi*4.6e9

#STIX coefficients
def getSDP(n, B):
    w_ce = -e*B/m_e
    w_cD = e*B/m_D
    w_pe2 = e**2*n/(m_e*eps_0)
    w_pD2 = e**2*n/(m_D*eps_0)

    S = 1 - (w_pe2/(w**2 - w_ce**2)) - (w_pD2/(w**2 - w_cD**2))
    D = (w_ce/w)*(w_pe2/(w**2 - w_ce**2)) + (w_cD/w)*(w_pD2/(w**2 - w_cD**2))
    P = 1 - w_pe2/w**2 - w_pD2/w**2
    return S, D, P


#determines which branch the waves is on by comparing the fast and slow solutions for N_perp^2 assuming a cold plasma
#to the N_perp^2 predicted by GENRAY. This therefore only works for cold plasma in GENRAY
#then plots either the ray trace (good for one or more rays) or the N_|| evolution (only readable for a single ray)
#plotEvolution -> N|| evo
#plotRayTrace -> raytrace
def getWhichRootPlots(plotEvolution = False, plotRayTrace = False):
    nparas = np.copy(genray_nc.variables["wnpar"]) #n_|| of the ray at each point along the ray trace
    nperps = np.copy(genray_nc.variables["wnper"]) #n_perp of the ray at each point along the ray trace
    radialVariable = (np.copy(genray_nc.variables["spsi"])) #rho_pol of the ray at each point along the ray trace
    wr  = genray_nc.variables["wr"][:]/100 #major radius of the ray at each point along the trace
    wr = np.ma.getdata(wr)
    wz  = genray_nc.variables["wz"][:]/100 #height of the ray at each point along the trace
    wz = np.ma.getdata(wz)
    ws  = genray_nc.variables["ws"][:]/100 #poloidal distance along ray
    delpwr= np.copy(cqlrf_nc.variables["delpwr"]) #power in the ray at each point

    #dielectric tensor components along the ray
    cweps11 = genray_nc.variables["cweps11"][:]
    cweps11 = cweps11[0] + 1j*cweps11[1]
    cweps12 = genray_nc.variables["cweps12"][:]
    cweps12 = cweps12[0] + 1j*cweps12[1]
    cweps21 = genray_nc.variables["cweps21"][:]
    cweps21 = cweps21[0] + 1j*cweps21[1]
    cweps33 = genray_nc.variables["cweps33"][:]
    cweps33 = cweps33[0] + 1j*cweps33[1]

    figEvo,axEvo = None, None
    figTrace, axTrace = None, None

    if plotEvolution:
        figEvo,axEvo = plt.subplots()#figsize = (7,4))#(6.4,4.8))
        axEvo.plot([-100], [-100], label = 'Fast Wave', lw = 2, color = 'r')
        axEvo.plot([-100], [-100], label = 'Slow Wave', lw = 2, color = 'b')

        axEvo.set_xlabel(r"$\rho_{pol}$")#; axEvoes[1].set_xlabel(r"$\rho_{pol}$")
        axEvo.set_ylabel(r'$N_{||}$')
        axEvo.set_ylabel(r'$N_{\perp}$')
        #axEvo.set_xlim([0,1.05])

    if plotRayTrace:
        figTrace,axTrace = plt.subplots(figsize = (5.25,7.1))
        plt.subplots_adjust(left=0.22,bottom = .1)
        axTrace.set_ylabel("z (m)")
        axTrace.set_xlabel("R (m)")

        #fake curves for legend
        axTrace.plot([-100], [-100], label = 'Fast Wave', lw = 2, color = 'r')
        axTrace.plot([-100], [-100], label = 'Slow Wave', lw = 2, color = 'b')

        xlim = gfileDict["xlim"] #R points of the wall
        ylim = gfileDict["ylim"] #Z points of the wall
        rbbbs = gfileDict["rbbbs"] #R points of the LCFS
        zbbbs = gfileDict["zbbbs"] # Z points of the LCFS

        axTrace.plot(xlim, ylim, 'grey', lw = 2)#plot wall
        axTrace.plot(rbbbs, zbbbs, 'k', lw = 1.5)#plot LCFS

        #ax.set_title(f"Plotting Rays until {(maxDelPwrPlot) * 100} %\n ray power deposition")
        axTrace.set_aspect('equal')
        
        axTrace.set_ylim(min(ylim)-.05, max(ylim)+.05)
        axTrace.set_xlim(min(xlim)-.05, max(xlim)+.05)
        axTrace.legend(loc = 'upper right')

        helper.drawFluxSurfaces(axTrace)

    norm = plt.Normalize(0, 1)
    minRatioToPlot = .05#.1 #when the ratio of ray power to ray starting power is below this number, the trace ends
    for i in range(len(nparas)):
        npara = nparas[i]
        nperp = nperps[i]
        
        N = np.sqrt(npara**2 + nperp**2)

        delpwrRatios = delpwr[i]/np.max(delpwr[i])
        mostPowerDep = helper.findNearestIndex(minRatioToPlot, delpwrRatios)

        S, D, P = np.real(cweps11[i]), np.imag(cweps21[i]), np.real(cweps33[i])
        
        #"""
        A = S
        B = -((S - npara**2)*(P + S) - D**2)
        C = P*((S - npara**2)**2 - D**2)
        Delta = B**2-4*A*C

        #this is the full, cold plasma accessibility condition
        #the equation that is much more often quoted is an approximation and gives the wrong answer
        N_acc = np.sqrt((-D**2*(P+S) + 2*np.sqrt(D**2*P*S*(D**2-(P-S)**2))+S*(P-S)**2)/((P-S)**2))

        #fast wave solution for N_perp^2
        fastNperp2 = (-B - np.sqrt(Delta))/(2*A)
        #slow wave solution for N_perp^2
        slowNperp2 = (-B + np.sqrt(Delta))/(2*A)
        
        fastNperp = np.sqrt(fastNperp2)
        fastNperp[fastNperp2 < 0] = np.inf
        slowNperp = np.sqrt(slowNperp2)
        slowNperp[slowNperp2 < 0] = np.inf
        
        distToFast = np.abs(np.abs(nperp)- fastNperp)
        distToSlow = np.abs(np.abs(nperp) - slowNperp)
        
        slowOrFast = np.zeros(len(npara))#zero where the wave is slow
        #if the nperp is closer to that of the fast branch, flip the relevant value to 1
        slowOrFast[distToFast < distToSlow] = 1

        if plotRayTrace:
            points = np.array([wr[i][:mostPowerDep], wz[i][:mostPowerDep]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            # Create a continuous norm to map from data points to colors
            lc = LineCollection(segments, norm = norm,cmap=plt.cm.bwr)
            # Set the values used for colormapping
            lc.set_array(slowOrFast[:mostPowerDep])
            lc.set_linewidth(2)
            axTrace.add_collection(lc)
            figTrace.tight_layout()

        if plotEvolution:
            #choose which ray to plot
            #There's not a great way to know which ray number is useful to look at besides guess and check
            if i is not 19:
                continue

            rho_pols = radialVariable[i]
            axEvo.plot(rho_pols[:mostPowerDep], np.sign(npara[0])*N_acc[:mostPowerDep], label = r'$N_{||, acc}$', lw = 4, color ='tab:green', linestyle = 'dotted')
            points = np.array([rho_pols[:mostPowerDep],npara[:mostPowerDep]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            # Create a continuous norm to map from data points to colors
            lc = LineCollection(segments, norm = norm,cmap=plt.cm.bwr)
            # Set the values used for colormapping
            lc.set_array(slowOrFast[:mostPowerDep])
            lc.set_linewidth(2.5)
            axEvo.add_collection(lc)
            axEvo.legend(loc = 'lower left')

            if nparas[0,0] > 0:
                axEvo.set_ylim([0,5])
            else:
                axEvo.set_ylim([-2.7,-2.4])
            axEvo.set_xlim([.81,.88])
            axEvo.set_yticks([-2.4,-2.5,-2.6,-2.7])

            axEvo.text(.866,-2.495, 'Propagation', rotation = -20, fontsize = 14, ha='center', va='center')
            axEvo.text(.871,-2.419, 'HFS', rotation = 0, fontsize = 14, ha='center', va='center')

            arrowStarts = [(.8565, -2.472), (.8312,-2.5215), (.859,-2.6068), (.866,-2.429)]
            arrowEnds= [(.8458,-2.452), (.842,-2.5512), (.856,-2.667), (.876,-2.429)]

            #arrow code for the paper figure
            linestyles = ['-','-','-','-']
            arrowStyles = ['simple', 'simple', 'simple', '->']
            widths = [3, 3, 3, 2]
            length = .0075
            print(f'before arrows: {radialVariable[i][mostPowerDep]}')
            for l in range(len(arrowStarts)):
                axEvo.annotate(
                        "",
                        xy=(arrowEnds[l][0], arrowEnds[l][1]),
                        xytext=(arrowStarts[l][0], arrowStarts[l][1]),
                        arrowprops=dict( linewidth = widths[l], 
                                        linestyle = linestyles[l], 
                                        color = 'k',
                                        arrowstyle = arrowStyles[l],
                                        joinstyle='miter',   # sharp corners
                                        capstyle='butt'  
                                        )
                        )
            figEvo.tight_layout()

    plt.show()

getWhichRootPlots(plotEvolution=False, plotRayTrace=True)
