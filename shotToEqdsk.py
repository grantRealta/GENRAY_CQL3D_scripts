###
# A dictionary for each machine mapping the shot to the appropriate gfile
###

Hammir_eqdskDict = {
    'baseCase' : 'eqdsk',
    'updated' : 'Hammir.eqdsk'
                    }

WHAM_eqdskDict = {
    'ECHtest' : 'WHAM.HTS_eqdsk',
    #'HTS' : 'WHAM.HTS_eqdsk',
    'HTS' : 'WHAM.HTS_eqdsk'
                  }



#shotNum is a string
def getEqdskName(shotNum, machine):
    if machine == 'Hammir':
        return Hammir_eqdskDict[shotNum]   
    if machine == 'WHAM':
        return WHAM_eqdskDict[shotNum]
