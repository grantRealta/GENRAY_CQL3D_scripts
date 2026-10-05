machine = 'WHAM'

if machine == 'Hammir':
    shot = 'updated'
    timeString = ''

    intermediaryDir = ''

    topmostShotDir = f'/global/homes/g/grutherf/scratch/Hammir_batch/{machine}_{shot}{timeString}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{timeString}'

if machine == 'WHAM':
    shot = 'HTS'
    timeString = ''

    intermediaryDir = '35GHz/'

    topmostShotDir = f'/global/homes/g/grutherf/scratch/WHAM_batch/{machine}_{shot}{timeString}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{timeString}_0.67zst_0.15rst_-46betast_1e19_40eV'
    """
    intermediaryDir = ''

    topmostShotDir = f'/global/homes/g/grutherf/scratch/WHAM_batch/{machine}_{shot}{timeString}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{timeString}_110GHzTest_clean'
    """

    """
    intermediaryDir = '110GHz_X2_1e19_40eV/'

    topmostShotDir = f'/global/homes/g/grutherf/scratch/WHAM_batch/{machine}_{shot}{timeString}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{timeString}_0.6zst_0.075rst_57betast_1e19_40eV'
    """


#returns the directory where the shot of interest is stored
def getTargetDir():
    return targetDir.strip()

#returns the machine where the shot is from (DIIID, WEST, etc)
def getMachine():
    print(f'machine: {topmostShotDir.split('/')[-1].split('_')[-2]}')
    return topmostShotDir.split('/')[-1].split('_')[-2]

def getShotNum():
    print(f'shotNum: {topmostShotDir.split('_')[-1]}')
    return topmostShotDir.split('_')[-1]

def getTopmostShotDir():
    return topmostShotDir