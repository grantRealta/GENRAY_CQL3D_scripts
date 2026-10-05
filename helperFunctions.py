###
# Series of helper functions that are used is a variety of other scripts
###


import numpy as np
import getGfileDict
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline
import matplotlib.pyplot as plt
import getTargetInfo
import getInputFileDictionary
import netCDF4
import shotToEqdsk
import f90nml

m_e = 9.109e-31
m_D = 3.343e-27
q_e = -1.602e-19
q_D = 1.602e-19
eps_0 = 8.854e-12

#this assumes that 
def getSPA_GEN(targetDir):
    if targetDir == None:
            targetDir = getTargetInfo.getTargetDir()

    genrayPath = f'{targetDir}/genray.in'
    nml = f90nml.read(genrayPath)
    assert nml['wave']['ireflm'] == 1

    gen_nc = netCDF4.Dataset(f'{targetDir}/genray_ec.nc','r')

    startingDelpwr = 0
    endingDelpwr = 0

    delpwr= gen_nc.variables["delpwr"][:] #power in the ray at each point

    for ray in range(len(delpwr)):
        startingDelpwr += delpwr[ray][0]
        #it seems like rays have their power set to zero when they are killed
        #so only look where the power is nonzero
        #this should be a rounding error for rays that deposit all their power
        endingDelpwr += (delpwr[ray][delpwr[ray]>0])[-1]

    return 1 - endingDelpwr/startingDelpwr
    

#returns the SPA of the forward lobe
def getSPA(targetDir):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()

    genray_in = None
    genray_in = getInputFileDictionary.getInputFileDictionary('genray_LH', targetDir=targetDir)
    cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')

    ngrill = genray_in['grill']['ngrill']
    lobes = np.arange(1,ngrill+1).astype(int)

    onePassDelpwrs= np.zeros(len(lobes))
    raysInLobes = np.zeros(len(lobes))
    startingDelpwrs = np.zeros(len(lobes))
    try:
        delpwr= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point
        radialVariable = (np.copy(cqlrf_nc.variables["spsi"]))
        nparas = np.ma.getdata(cqlrf_nc.variables['wnpar'][:])
        for i in range(len(lobes)):
            lobe = lobes[i]

            if lobe > genray_in['grill']['ngrill']:
                continue
                
            for ray in range(len(delpwr)):
            
                if genray_in['grill'][f'anmax({lobe})'] >= nparas[ray][0] >= genray_in['grill'][f'anmin({lobe})']:
                    firstBounceIndex = findBounceIndex(radialVariable[ray],bounceToFind = 1)
                    
                    startingDelpwrs[i] += delpwr[ray][0]
                    onePassDelpwrs[i] += delpwr[ray][firstBounceIndex]
                    raysInLobes[i] +=1
        avgSPA = 1-(onePassDelpwrs/startingDelpwrs)
        avgSPA[avgSPA < 0] = np.nan
        avgSPA[avgSPA > 1] = np.nan

        return avgSPA, onePassDelpwrs, startingDelpwrs
    
    except:
        import traceback
        traceback.print_exc()
        print(f'something went wrong')
        return -1

#returns the index of the array whose element is closest to value
def findNearestIndex(value, array):
    idx = (np.abs(array - value)).argmin()

    return idx

#returns the index at which the input ray bounces in the SOL
def findBounceIndex(radius, bounceToFind = 1):
    checkpointValue = 1
    pastCheckpoint = False

    numBounces = 0
    for i in range(len(radius)):
        if radius[i] <= checkpointValue:
            pastCheckpoint = True
        if radius[i] > checkpointValue and pastCheckpoint:
            numBounces += 1
            pastCheckpoint = False
            if numBounces == bounceToFind:
                return i
    return i

#returns the location of where the driven current is maximum
#this may or may not line up with the where the current density is highest
def getLocationOfCurrentPeak(cql_nc):
    import matplotlib.pyplot as plt
    curr = cql_nc.variables["curr"][-1,:]
    dArea = cql_nc.variables["darea"][:]
    paraCurrentProf = (curr * dArea)
    rya = cql_nc.variables["rya"][:]
    return rya[np.argmax(np.abs(paraCurrentProf))]

#returns the location of where the driven current density is maximum
#this may or may not line up with the where the current is highest
def getLocationOfCurrentDensityPeak(cql_nc):
    import matplotlib.pyplot as plt
    curr = cql_nc.variables["curr"][-1,:]
    dArea = cql_nc.variables["darea"][:]
    paraCurrentProf = (curr * dArea)
    rya = cql_nc.variables["rya"][:]
    return rya[np.argmax(np.abs(curr))]

#draw poloidal flux surfaces
def drawFluxSurfaces(ax, targetDir = None, rhosToPlot = [.1, .5, 1], 
                     colors = 'k', zBounds = None, limPath = None, bothSides = True):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    
    gfileDict = getGfileDict.getGfileDict(targetDir=targetDir)
    r = gfileDict["rgrid"]
    z = gfileDict["zgrid"]
    psirz = gfileDict["psirz"]
    
    psi_mag_axis = gfileDict["ssimag"]
    psi_boundary = gfileDict["ssibdry"]
    
    if zBounds == None:
        zBounds = [np.min(z), np.max(z)]
    psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)

    rInterp = np.linspace(np.min(r), np.max(r), 201)
    zInterp = np.linspace(zBounds[0], zBounds[1], 200)

    interp = RectBivariateSpline(r, z, psirzNorm.T, kx=3, ky=3)
    psirzNormInterp = interp(rInterp, zInterp)

    if limPath is not None:
        Rs, Zs = np.meshgrid(rInterp, zInterp)
        points = np.vstack((Rs.ravel(), Zs.ravel())).T
        mask = limPath.contains_points(points).reshape(Rs.shape)
        psirzNormInterp[np.logical_not(mask)] = np.nan     

    ax.contour(zInterp, rInterp, psirzNormInterp, np.square(rhosToPlot), colors= colors, linewidths=2.25)

    if bothSides:
        ax.contour(zInterp, -rInterp, psirzNormInterp, np.square(rhosToPlot), colors= colors, linewidths=2.25)

#returns rya and the absolute value of the safety factor profile
def getCQLq(targetDir = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
    rya = cql_nc.variables["rya"][:]
    q_prof = cql_nc.variables["qsafety"][:]

    return rya, np.abs(q_prof)

#returns the electron density profile in cql3d
def getCQLne(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    enescal = 1
    try:
        enescal = cqlInputDict['setup']['enescal']
    except:
        print(f'enescale not in input dict')
    n_e = cqlInputDict['setup']['enein(1,1)']*1e6*enescal
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        n_e = interp1d(rhos, n_e)(rho_pol)
        rhos = rho_pol

    return rhos, n_e

#returns the CQL3D Te profile
def getCQLTe(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    tescal = 1
    try:
        tescal = cqlInputDict['setup']['tescal']
    except:
        print(f'tescal not in input dict')
    T_e = cqlInputDict['setup']['tein']*tescal
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        T_e = interp1d(rhos, T_e)(rho_pol)
        rhos = rho_pol

    return rhos, T_e

#returns the deuterium density profile in cql3d
def getCQLnD(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    enescal = 1
    try:
        enescal = cqlInputDict['setup']['enescal']
    except:
        print(f'enescale not in input dict')
    n_D = cqlInputDict['setup']['enein(1,2)']*1e6*enescal
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        n_D = interp1d(rhos, n_D)(rho_pol)
        rhos = rho_pol

    return rhos, n_D

#returns the electron density profile in cql3d
def getCQLnC(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    enescal = 1
    try:
        enescal = cqlInputDict['setup']['enescal']
    except:
        print(f'enescale not in input dict')
    n_C = cqlInputDict['setup']['enein(1,3)']*1e6*enescal
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        n_C = interp1d(rhos, n_C)(rho_pol)
        rhos = rho_pol

    return rhos, n_C

#returns the CQL3D Ti profile
def getCQLTD(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    tiscal = 1
    try:
        tiscal = cqlInputDict['setup']['tiscal']
    except:
        print(f'tiscal not in input dict')
    T_i = cqlInputDict['setup']['tiin']*tiscal
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        T_i = interp1d(rhos, T_i)(rho_pol)
        rhos = rho_pol

    return rhos, T_i

#returns the CQL3D Zeff profile
def getCQLZeff(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    Zeff = cqlInputDict['setup']['zeffin']
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        Zeff = interp1d(rhos, Zeff)(rho_pol)
        rhos = rho_pol

    return rhos, Zeff

#convert the input rho_pol coordinate to major radius at the midplane
#you can choose whether you take the HFS or LFS midplane values
def convertRhopolToRmidplane(rhos, targetDir = None, side = 'LFS'):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    gfileDict = getGfileDict.getGfileDict(targetDir=targetDir)

    r = gfileDict["rgrid"]
    z = gfileDict["zgrid"]
    psirz = gfileDict["psirz"]
    
    psi_mag_axis = gfileDict["ssimag"]
    psi_boundary = gfileDict["ssibdry"]
    
    psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)

    rInterp = np.linspace(np.min(r), np.max(r), 200)
    psirzNormFunc = interp2d(r,z, psirzNorm, kind = 'cubic')

    zmaxis = gfileDict["zmaxis"]
    Rmaxis = gfileDict['rmaxis']

    if side == 'LFS':
        relevantR_interp = rInterp[rInterp >= Rmaxis]
    elif side == 'HFS':
        relevantR_interp = rInterp[rInterp <= Rmaxis]
    else:
        relevantR_interp = rInterp
        print(f'returning both HFS and LFS values')

    midplanePsi = psirzNormFunc(relevantR_interp, zmaxis)
    
    Rs_fromRho_psi = interp1d(midplanePsi, relevantR_interp, fill_value = 'extrapolate', bounds_error = False)(rhos**2)
    return Rs_fromRho_psi

#convert the input major radius coordinate to rho_pol
def convertRtoRhopol(Rs, targetDir = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    gfileDict = getGfileDict.getGfileDict(targetDir=targetDir)
    r = gfileDict["rgrid"]
    z = gfileDict["zgrid"]
    psirz = gfileDict["psirz"]
    
    psi_mag_axis = gfileDict["ssimag"]
    psi_boundary = gfileDict["ssibdry"]
    
    psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)

    rInterp = np.linspace(np.min(r), np.max(r), 200)
    psirzNormFunc = interp2d(r,z, psirzNorm, kind = 'cubic')

    zmaxis = gfileDict["zmaxis"]
    Rmaxis = gfileDict['rmaxis']

    midplanePsi = psirzNormFunc(rInterp, zmaxis)
    
    rhopols_from_Rs = np.sqrt(interp1d(rInterp, midplanePsi, fill_value = 'extrapolate', bounds_error = False)(Rs))
    return rhopols_from_Rs
    
def get_wpe2(ne):
    return ne*q_e**2/(m_e*eps_0)

def get_wpD2(ne):
    return ne*q_D**2/(m_D*eps_0)

def get_Wce(B):
    return q_e*B/m_e

def get_WcD(B):
    return q_D*B/m_D