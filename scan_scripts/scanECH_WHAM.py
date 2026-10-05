import numpy as np
import os, sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import setupInputFiles
import getGfileDict


rho_gen = np.linspace(0, 1, 51)
n_0 = 1e19; n_1 = 1e18; r_1 = 3.2; r_2 = .4
den_prof = (n_0-n_1)*(1-rho_gen**r_1)**r_2+n_1
den_prof = [float(f"{elem:.2e}") for elem in den_prof]

n_e = np.copy(den_prof); n_i=np.copy(den_prof)
dentab = [x for tup in zip(n_e, n_i) for x in tup]

T_0 = .04; T_1 = .01; t_1 = 3.2; t_2 = .4
tem_prof = (T_0-T_1)*(1-rho_gen**t_1)**t_2+T_1
tem_prof = [float(f"{elem:.2e}") for elem in tem_prof]

T_e = np.copy(tem_prof); T_i=np.copy(tem_prof)
temtab = [x for tup in zip(T_e, T_i) for x in tup]


rho_cql = np.linspace(0, 1, 101)
n_0 = 1e19; n_1 = 1e18; r_1 = 3.2; r_2 = .4
den_prof_cql = (n_0-n_1)*(1-rho_cql**r_1)**r_2+n_1

den_prof_gen = interp1d(rho_cql, den_prof_cql)(rho_gen)

dentab = [x for tup in zip(den_prof_gen, den_prof_gen) for x in tup]

T_0 = .04; T_1 = .01; t_1 = 3.2; t_2 = .4

Tscale = 1

T_prof_cql = (T_0-T_1)*(1-rho_cql**t_1)**t_2+T_1
T_prof_cql *= Tscale

T_prof_gen = interp1d(rho_cql, T_prof_cql)(rho_gen)

temtab = [x for tup in zip(T_prof_gen, T_prof_gen) for x in tup]

genEdits = {
    'wave':{
        'frqncy':35e+09,
        'ioxm': -1
    },
    'eccone':{
        'zst' : .75,
        'rst' : .7
    },
    'dentab':{
        'prof':dentab
    },
    'temtab':{
        'prof':temtab
            }
    }
den_prof_cql = (den_prof_cql/1e6).tolist()
cqlEdits = {
    'setup':{
        'enein':[den_prof_cql, den_prof_cql, den_prof_cql],
        'tein':T_prof_cql,
        'tiin' :T_prof_cql,
        'ryain':rho_cql,
    },
    
    }
intermediateDir = '35GHz/'
targetStem = f'/global/homes/g/grutherf/codes/IPS_scripts/WHAM_batch/WHAM_HTS/'
gfileDict = getGfileDict.getGfileDict(targetDir = targetStem)
targetStem = f'{targetStem}{intermediateDir}/WHAM_HTS'

rlim = gfileDict["xlim"] #R points of the wall
zlim = gfileDict["ylim"] #Z points of the wall

rlimFunc = interp1d(zlim, rlim)

rst = .15#0.08#.066
zsts = [.67]#[.6,.65,.7]#,.7,.8]
for i in range(len(zsts)):
    zst = zsts[i]

    print(f'zst: {zst}, rst: {rst}')

    center_betast = -np.degrees(np.atan((zst-.56)/(rst-0)))
    closeEdge_betast = -np.degrees(np.atan((zst-.56)/(rst - .045)))

    print(f'zst: {zst}, rst: {rst}, center betast: {center_betast}, close betast: {closeEdge_betast}')

    betasts = np.linspace(center_betast, closeEdge_betast,10)
    #betasts = [-70,-80,-90]
    print(betasts)

    for j in range(len(betasts)):
        betast = betasts[j]
        betast = np.round(betast)
        print(betast)
    #"""
        genEdits['eccone']['zst']=zst
        genEdits['eccone']['rst']=rst
        genEdits['eccone']['betast']=betast


        targetDir = f'{targetStem}_{zst}zst_{rst}rst_{betast:.0f}betast_1e19_40eV'

        inputFileHelper = setupInputFiles.InputFileHelper(targetDir = targetDir, genray_edits=genEdits,cql_edits = cqlEdits)
        inputFileHelper.copyInputFileTemplates()
        inputFileHelper.setupGENRAY()
        inputFileHelper.setupCQL()

        os.system(f'cd "{targetDir}" && bash runGEN.sh "{targetDir}"')
    #"""