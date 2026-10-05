#GENRAY makes a lot of files. This script removes many of the ones that aren't necessary for data analysis
#I suggest running this regularly but be sure not to run it while the codes are also running


import numpy as np
import os, sys
import netCDF4
import matplotlib.pyplot as plt


def garbageCollection():
    rootDir = f'/home/grantr/symlinks/genray_batch/DIIID_shots/'
    folderList = fast_scandir(rootDir)
    
    fig, ax = plt.subplots()

    mostCountsDir = ''
    mostCounts = -1

    for folder in folderList:
        print(folder)
        try:
            os.system(f'rm {folder}/*.bin')
            os.system(f'rm {folder}/*.doc')
            os.system(f'rm {folder}/*.sap')
            os.system(f'rm {folder}/*.dat')
            os.system(f'rm {folder}/con1')
            os.system(f'rm {folder}/genray.txt')
            os.system(f'rm {folder}/genray_one_ray_point.nc')
        except:
            pass

def fast_scandir(dirname):
    subfolders= [f.path for f in os.scandir(dirname) if f.is_dir()]
    for dirname in list(subfolders):
        subfolders.extend(fast_scandir(dirname))
    return subfolders  


garbageCollection()
