# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 11:33:22 2018

Uses the edf2asc from SR research function "EyeLinkDevKit_Windows_1.11.5.exe"

@author: lindadevoogd
"""

import sys
import os
import glob
import subprocess
import shutil

# SETTINGS
#------------------------------------------------------------------------------

#paths
homepath="D:\\"

projectpath=os.path.join(homepath,'NYU_RS_LC')
rawdatapath=os.path.join(projectpath,'rawdata')
datapath=os.path.join(projectpath,'data')
statspath=os.path.join(projectpath,'stats')
sys.path.append(os.path.join(projectpath,'scripts'))

def edf2asc_py(edffile):
    
    # Download the edf2asc tool from the website:
    # https://www.sr-research.com/downloads/
    #
    # import subprocess
    
    # Covert the EDF to ASC
    cmdCommand = "edf2asc " + edf_file   #specify your cmd command
    process = subprocess.Popen(cmdCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
for filename in glob.glob(os.path.join(rawdatapath,"MRI_FCWML*")):
    
    for days in "logfiles-ses-day1","logfiles-ses-day2":
        
        edf_files=glob.glob(os.path.join(filename,days,"*.EDF"))
        
        for edf_file in edf_files:

            # Covert the EDF to ASC
            edf2asc_py(edf_file)

            # Move .ASC file to data
            subjname = os.path.basename(filename)
            pupfilename = os.path.basename(edf_file)
            os.mkdir(os.path.join(datapath,subjname,days))
            despath=os.path.join(datapath,subjname,days,pupfilename[:-4] + '.asc')
            shutil.move(edf_file[:-4] +'.asc', despath)
    




