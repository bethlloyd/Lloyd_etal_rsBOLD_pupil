# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 13:37:21 2020

@author: linda
"""

# Import
from EyeLinkRead import EyeLinkRead # own model in scripts folder
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
import scipy.ndimage as ndimage
import re
import os
import glob
import numpy as np

# Resting state session contains 150 .nii images with a TR of 2 sec, the total 
# duration should be 300 sec. Each session has 5 dummy scans and this preceeds 
# the task script start as well as the "START_TIME" marker in the eye link data
# In the task script this is the sequence of events:
#
# 5 dummy scans = 10 sec
# start of eye link (this causes a delay of 375 ms...)
# sending of "START_TIME" marker to eye link
# wait 2000 ms (this was done for task purposes)
# logging of start_time in log file of the task
# RS duration of 300 sec
# logging of end_time in log file of the task > thus end_time - start_time is 300 sec
# 2000 ms of the goodbye text
# sending of "END_TIME" marker to eye link > this END_TIME - START_TIME should be 304 sec
# stop recording and close eye link

#
# > there seems to be ~500 ms more between START_TIME and END_TIME than it should be
#


def downsample2bins(timecoursedata,windowlength,nbins):
    
    """
    This function splits up timecourse data into average bins
    
    Input:
        -timecoursedata: your timecourse data
        -windowlength: length of the bin you want an average off [must be in the
        same frequency as the timecourse]
        -nbins: number of bins you want to average the data to
    Output:
        -list with the averages of each bin
    
    """
    
    # Check if lengths match up
    if len(timecoursedata)==windowlength*nbins:
    
        # Make an empty vector
        list_bins=[]
        
        # Window settings for the first bin
        start_pos=0
        end_pos=windowlength
        
        # Loop over the number of bins
        for c in range(0,nbins):
            
            if c==0:
                
                # For the firts bin take the average
                list_bins.append(np.average(timecoursedata[start_pos:end_pos]))
                
                
            else:
                
                # For all but first bin, use the updated window
                list_bins.append(np.average(timecoursedata[start_pos:end_pos]))
                
            # Update the window for the next round
            start_pos=end_pos+1
            end_pos=end_pos+windowlength
            
        return list_bins
    
    else:
        
        print('The length of your timecourse does not match the windowlength and nBins')



# Proces pupil function
#------------------------------------------------------------------------------


def a_pupilprocessing(filename,nScans,TR):

    # Load in the class
    eyelink_module=EyeLinkRead()
    
    # Load in the file
    eyelink_module.do_readfile(filename)
    
    # Check position START_TIME, should be 93
    starttimeposition=[cnt for cnt, line in enumerate(eyelink_module.rawdat) if "START_TIME" in line]
    #print("START_TIME position is " + str(starttimeposition))
    
    # Cut out desired time course
    scanduration=nScans*TR*eyelink_module.sF
    starttime=0
    endtime=scanduration
    eyelink_module.do_cutoutdata(starttime,endtime)
    
    # Get eye blinks for interpolation
    eyelink_module.get_eyeblinks()
    
    # Interpolate eye blinks
    eyelink_module.do_interpol()
    
    # Smooth [change to your liking or leave this step out]
    eyelink_module.smooth_int_pupdat=ndimage.filters.gaussian_filter(eyelink_module.int_pupdat,100)
    
    # Now average the time course to bins based on TR length
    #windowlength=TR*eyelink_module.sF
    pup_list=downsample2bins(eyelink_module.int_pupdat,TR*sF, nScans)
    c_pup_list = [val-np.average(pup_list) for val in pup_list]
    eyelink_module.z_pup_list = [val/np.std(c_pup_list) for val in c_pup_list]
    
    # Get derivitive of the pupil
    eyelink_module.z_pup_list_diff=list(np.diff(eyelink_module.z_pup_list))
    eyelink_module.z_pup_list_diff.insert(0,0)

    return eyelink_module


        


# Run the pupil processing
#------------------------------------------------------------------------------

# Settings
nScans=150
TR=2

# Path settings
homepath = "C:\\Users\\lloydb\\surfdrive\\ExperimentData\\NYU_RS_LC"
raw_datapath = os.path.join(homepath, 'rawdata')
datapath = os.path.join(homepath, 'data')

# specify day
day = ['day1', 'day2']

# Loop over subjects
for subj in range(len(os.listdir(datapath))):

    # Loop over days
    for d in range(len(day)):
        
                
        # Name current file 
        filename = glob.glob(os.path.join(datapath, os.listdir(datapath)[subj], 'logfiles-ses-' + day[d],  "MRI*.asc"))
        #print(filename)

        # Run the preprocessing
        try:
            
            # Do the preproceccing
            eyelink_module = a_pupilprocessing(filename[0],nScans,TR)
            
            # Print out average surface
            #print("Average surface area is " + str(np.average(eyelink_module.int_pupdat)))
            
            # Print out % interpolated data
            interpol_vec=[int(x==0) for x in eyelink_module.pupdat]
            prop_int=sum(interpol_vec)/len(eyelink_module.pupdat)
            print("% of interpolated data is " + str(prop_int*100) + " %")
            
            # Plot interpolated data over rawdata
            plt.figure(1)
            plt.title('raw data')
            plt.plot(eyelink_module.pupdat, 'b')
            plt.xlabel('pupil samples (250Hz)')
            plt.ylabel('pupil size (area units)')
            plt.title('raw data (blue), processed data (red)')
            plt.plot(eyelink_module.int_pupdat, 'r')
            plt.plot(eyelink_module.smooth_int_pupdat, 'g')
            
            # Plot TR binned data + derivative
            plt.figure(2)
            plt.title('fmri time-locked data')
            plt.plot(eyelink_module.z_pup_list)
            plt.xlabel('time series (TR)')
            plt.ylabel('pupil size (z-score)')
            plt.plot(eyelink_module.z_pup_list_diff)
        
        
            # Save the data
            head_tail = os.path.split(filename[0]) 
            savepath=head_tail[0]
            
            np.savetxt(os.path.join(head_tail[0],"pupil_dilation.txt"), eyelink_module.z_pup_list, delimiter=',')
            np.savetxt(os.path.join(head_tail[0],"pupil_derivative.txt"), eyelink_module.z_pup_list_diff, delimiter=',')
            
        except IndexError: # catch the error
            print("data missing!")
            
            
            from itertools import islice

            def downsample_to_proportion(rows, proportion=1):
                return list(islice(rows, 0, len(rows), int(1/proportion)))

            
            #check something with timing: downsample to 50Hz
            pupdat50hz = downsample_to_proportion(eyelink_module.pupdat, 0.2)


    







#make lists
file_list = []
file_length = []
length_difference = []


# # Read in the .ASC file
# # go through each subject 
# for subj in range(len(os.listdir(datapath))):
#     # go through both days
#     for d in range(len(day)):
        
#         # Load in the class
#         eyelink_module=EyeLinkRead()
        
#         # name current file 
#         file = glob.glob(os.path.join(datapath, os.listdir(datapath)[subj], 'logfiles-ses-' + day[d],  "MRI*.asc"))
        
#         # read in file
#         try:
#             eyelink_module.do_readfile(file[0])
#             # log datafile length
#             #file_length.append(str(eyelink_module.endtime-eyelink_module.starttime))
#             # log difference in length
#             #length_difference.append(int(eyelink_module.endtime-eyelink_module.starttime)-304000)
#             #file_list.append(file)
#         except IndexError: # catch the error
#             pass # pass will basically ignore it
#         # log length of data file 
       
        
#         del eyelink_module
        

#
#eyelink_module.do_readfile('MRI_FCWML001day1rest1.asc')
#print("data length should be 304000  but is " + str(eyelink_module.endtime-eyelink_module.starttime))
##day1 304620
##day2 304504

