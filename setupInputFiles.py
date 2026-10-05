#!/usr/bin/env python3

import f90nml
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from scipy.interpolate import interp1d
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getTargetInfo
import shotToEqdsk

IPS_scripts_location = '/global/homes/g/grutherf/codes/IPS_scripts'

class InputFileHelper:
    """
    
    """
    def __init__(self, targetDir = None,  
                 genray_edits = None,
                 cql_edits = None
                ):

        if targetDir == None:
            self.targetDir = getTargetInfo.getTargetDir()
        else:
            self.targetDir = targetDir

        targetSplit = targetDir.split('/')
        idx = next(p for p, s in enumerate(targetSplit) if "batch" in s)
        self.parentShotDir = targetSplit[idx + 1]
        self.shot = self.parentShotDir.split('_')[1]
        self.machine = self.parentShotDir.split('_')[0]

        self.eqdskName = shotToEqdsk.getEqdskName(self.shot, self.machine)

        self.genray_edits = genray_edits
        self.cql_edits = cql_edits

        #self.cql_edits['eqsetup']={'eqsource':self.eqdskName}
        #self.genray_edits['tokamak']={'eqdskin':self.eqdskName}

    def copyInputFileTemplates(self):
        
        #if the target dir doesn't exist, either make one or break
        if not os.path.exists(self.targetDir):
            print(f'making dir {self.targetDir}')
            os.system(f'mkdir -p {self.targetDir}')

        print(f'copying in or overwriting input file templates')

        os.system(f'cp {IPS_scripts_location}/templates/{self.machine}_templates/cqlinput {self.targetDir}/cqlinput')
        os.system(f'cp {IPS_scripts_location}/templates/{self.machine}_templates/genray.in {self.targetDir}/genray.in')
        os.system(f'cp /global/homes/g/grutherf/codes/IPS_scripts/{self.machine}_batch/{self.parentShotDir}/{self.eqdskName} {self.targetDir}/eqdsk')

        os.system(f'cp {IPS_scripts_location}/runGENThenCQL.sh {self.targetDir}')
        os.system(f'cp {IPS_scripts_location}/runGEN.sh {self.targetDir}')
        os.system(f'cp {IPS_scripts_location}/runCQL.sh {self.targetDir}')
        os.system(f'cp {IPS_scripts_location}/cql.slurm {self.targetDir}')
        os.system(f'cp {IPS_scripts_location}/genr.slurm {self.targetDir}')
        os.system(f'cd {self.targetDir} && chmod +x *.sh')


    def setupGENRAY(self):
        genrayPath = f'{self.targetDir}/genray.in'

        nml = f90nml.read(genrayPath)

        for group in self.genray_edits.keys():
            for variable in self.genray_edits[group].keys():
                nml[group][variable] = self.genray_edits[group][variable]

        nml.write(genrayPath, force=True)
    
        print(f"Wrote GENRAY")

    def setupCQL(self):
        cqlPath = f'{self.targetDir}/cqlinput'

        nml = f90nml.read(cqlPath)
        for group in self.cql_edits.keys():
            for variable in self.cql_edits[group].keys():
                nml[group][variable] = self.cql_edits[group][variable]

        nml.write(cqlPath, force=True)
    
        print(f"Wrote CQL")


def main():
    targetDir = f'/global/homes/g/grutherf/codes/IPS_scripts/WHAM_batch/WHAM_HTS/WHAM_HTS_110GHzTest'

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
    T_prof_cql = (T_0-T_1)*(1-rho_cql**t_1)**t_2+T_1

    T_prof_gen = interp1d(rho_cql, T_prof_cql)(rho_gen)

    temtab = [x for tup in zip(T_prof_gen, T_prof_gen) for x in tup]


    genEdit = {
        'wave':{
            'frqncy':110e+09,
        },
        'eccone':{
            'zst' : .874,
            'rst' : .212
        },
        'dentab':{
            'prof':dentab
        },
        'temtab':{
            'prof':temtab
                }
        }
    print(f'den_prof_cql: {den_prof_cql/1e6}')
    den_prof_cql = (den_prof_cql/1e6).tolist()
    cqlEdits = {
        'setup':{
            'enein':[den_prof_cql, den_prof_cql, den_prof_cql],
            'tein':T_prof_cql,
            'tiin' :T_prof_cql,
            'ryain':rho_cql,
        },
        
        }
    inputFileHelper = InputFileHelper(targetDir = targetDir, genray_edits=genEdit,cql_edits = cqlEdits)
    inputFileHelper.copyInputFileTemplates()
    inputFileHelper.setupGENRAY()
    inputFileHelper.setupCQL()


if __name__ == "__main__":
    main()