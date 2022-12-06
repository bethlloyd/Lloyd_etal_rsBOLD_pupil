# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 12:38:50 2020

@author: lloydb
"""

# Import
#from EyeLinkRead import EyeLinkRead # own model in scripts folder
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
import scipy.ndimage as ndimage
from scipy import stats
import os
import glob
import numpy as np
from itertools import islice
import pandas as pd 
from scipy.stats import gamma
from sklearn.preprocessing import MinMaxScaler
def downsample_to_proportion(rows, proportion=1):
    return list(islice(rows, 0, len(rows), int(1/proportion)))
from scipy import signal
import stimfuncs
import statistics
#import nideconv as nd
#import pystan
#from nideconv import GroupResponseFitter
import seaborn as sns
# Settings
nScans=150
TR=2
sF=50 #output from by PupCor GUI

# averag signal settings: 
baseline_dur=1
epoch_dur = 10


# Functions
def downsample2bins(timecoursedata,windowlength,nbins,correct_sd):
    
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
            
            #if c==0:
                
            # For the firts bin take the average
            if correct_sd==1:
            
                # Remove samples +/1 3SD outside mean calculation [see Murphy et. al. 2014]
                pup_list = [abs(number) for number in timecoursedata[start_pos:end_pos]]
                temp_3sd=(np.std(pup_list)*3)
                temp_min=np.average(pup_list)-temp_3sd
                temp_plus=np.average(pup_list)+temp_3sd
                temp_timecourse=pup_list
                
                #check if the sample sits within min and max threshold 
                bool1=[val<temp_min for val in temp_timecourse]
                bool2=[val>temp_plus for val in temp_timecourse]
                bool3=np.logical_or(bool1, bool2)
                
                prop_invalid_samples=sum(bool3)/len(bool3)
                
                if prop_invalid_samples > 0.40:
                    print('prop invalid samples in epoch = ', prop_invalid_samples)
                # take the new average of the epoch and save to list 
                list_bins.append(np.average([val for c,val in enumerate(temp_timecourse) if bool3[c]==False]))
                

            elif correct_sd==0:

                list_bins.append(np.average(timecoursedata[start_pos:end_pos]))
                
            # Update the window for the next round
            start_pos=end_pos+1
            end_pos=end_pos+windowlength
            
        return list_bins
    
    else:
        
        print('The length of your timecourse does not match the windowlength and nBins')

def get_invalid_samples(filename,sF):
    
    f = open(filename, 'r')
    rawdat=f.readlines()
    
    # Get the sample frequency
    rcd_line=[line for cnt, line in enumerate(rawdat) if "RECCFG" in line]
    sampleline=rcd_line[0].split()
    eyelink_sF=int(sampleline[4])
    print('Your sample freq is ' + str(eyelink_sF) + ' Hz and will be downsamled to 50 Hz')

    # remove lines [script crashes when these lines are in there]
    for rmstr in ["EFIX","ESACC","SFIX","SSACC","SBLINK","EBLINK","END"]:
        rawdat = [line for line in rawdat if not rmstr in line]
    
    # take the pupil diameter
    dat=[ [], [], [], [] ]
    rawevt=[]
    
    #get recording line so everything  before will be removed
    rcd_line=[cnt for cnt, line in enumerate(rawdat) if "SAMPLES" in line]
    
    # remove those lines from rawdata
    rawdat = rawdat[rcd_line[-1]+1:]
    
    # get events and make vector
    for c,line in enumerate(rawdat):
        #if c>rcd_line[-1]: # just throws away lines before calibration [but might be file specific]

        if "MSG" in line:
            rawevt.append(line)
        else:
            spt_line=line.split('\t')
            for cc,dd in enumerate(spt_line):
                if cc<4:
                    dat[cc].append(dd.strip('  '))
    
    #get pupil dilation
    pupdat = [int(float(x)) for x in dat[3]]

    #down sample to 50 HZ
    downsF=int(eyelink_sF/sF)
    pupdat=pupdat[0::downsF]
    
    prop=(len([val for val in pupdat if val==0])/len(pupdat))*100
    print("% of invalid samples is " + str(prop) + " %")  
    return prop



def get_point_processEvents(timcourse, threshold, sF):   
    # timecourse = signal (as list), threshold (no. of sd to set as threshold), sF = 50
    
    centered_TC = timcourse-np.average(timcourse)
    sd_pup = statistics.pstdev(centered_TC)
    mean_pup = statistics.mean(centered_TC)
    
    cutoff_up = mean_pup+(sd_pup*threshold)
    
     #check if the sample sits within min and max threshold 
    bool1=[val>cutoff_up for val in centered_TC]

    threshold_signal = []
    for count, x in enumerate(centered_TC):
        if bool1[count] == True: 
            threshold_signal.append(x)
        else:
            threshold_signal.append(sd_pup)
    
    # get onsets of psuedo events
    onsets = []
    for count, y in enumerate(threshold_signal):
        if count < len(threshold_signal)-1:
            if y == sd_pup:
                if threshold_signal[count+1] > sd_pup:
                    #get index 
                    onsets.append((count-sF)/sF) # convert to s
                    
    
    return onsets
    

# cut out sample of data 
def extract_ave_epoch(timecourse, onsets, baseline_dur, epoch_dur, sF):

    agg_epoch = make_empty_array(0, baseline_dur+epoch_dur, sF)
    # get start sample and end sample 
    epoch_start  = [(x-baseline_dur)*sF for x in BOLD_onsets]
    epoch_end = [(x+epoch_dur)*sF for x in BOLD_onsets]
    
    for index, (start, end) in enumerate(zip(epoch_start, epoch_end)):
        if end >30000:
            agg_epoch=np.vstack([agg_epoch, make_empty_array(0, baseline_dur+epoch_dur, sF)])
        else:
            agg_epoch=np.vstack([agg_epoch,  timecourse[start: end]])
    if len(agg_epoch) != len(onsets+1):
        print('oi! incorrect number of events happening somewhere')
    
   # epochs_mean = np.nanmean(agg_epoch[:],axis=0)
    return agg_epoch



def make_empty_array(time_start, time_end, sF):
    #empty_array = np.array(list(range(1,int((time_end-time_start)*sF)+1)))
    empty_array = np.empty((1,int((time_end--time_start)*sF)))
    empty_array[:] = np.NaN # make array list 
    empty_array=np.array(empty_array).ravel() # unravel
    
    return empty_array


    
from scipy.signal import butter, sosfilt, sosfreqz

def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfilt(sos, data)
        return y


#------------------------------------------------------------------------------
# Path settings
homepath = "D:\\NYU_RS_LC"
raw_datapath = os.path.join(homepath, 'data')
group_HRF_path = os.path.join(homepath, 'stats', 'group_HRFs')
# to add later:
ROI_HRF_path ="F:\\NYU_RS_LC\\stats\\rsHRF\\groupstats\\"
#------------------------------------------------------------------------------

# load in the SPM canonical HRF
SPM_canonicalHRF = os.path.join(group_HRF_path, 'SPM_canonicalHRF_plus.csv')
SPM_canonicalHRF = pd.read_csv(SPM_canonicalHRF, sep=';', header=None)

# load in the SPM canonical HRF - with default peak at 5s 
SPM_canonicalHRF_P5 = os.path.join(group_HRF_path, 'SPM_canonicalHRF_shift5.csv')
SPM_canonicalHRF_P5 = pd.read_csv(SPM_canonicalHRF_P5, sep=';', header=None)
SPM_canonicalHRF_P5 = SPM_canonicalHRF_P5.loc[:, 0]

# load in the SPM canonical HRF - with default peak at 4s (instead of 6s) 
SPM_canonicalHRF_P4 = os.path.join(group_HRF_path, 'SPM_canonicalHRF_shift4.csv')
SPM_canonicalHRF_P4 = pd.read_csv(SPM_canonicalHRF_P4, sep=';', header=None)
SPM_canonicalHRF_P4 = SPM_canonicalHRF_P4.loc[:, 0]

# load in the SPM canonical HRF - with default peak at 3s (instead of 6s) 
SPM_canonicalHRF_P3 = os.path.join(group_HRF_path, 'SPM_canonicalHRF_shift3.csv')
SPM_canonicalHRF_P3 = pd.read_csv(SPM_canonicalHRF_P3, sep=';', header=None)
SPM_canonicalHRF_P3 = SPM_canonicalHRF_P3.loc[:, 0]

# load in the SPM canonical HRF - with default peak at 2s (instead of 6s) 
SPM_canonicalHRF_P2 = os.path.join(group_HRF_path, 'SPM_canonicalHRF_shift2.csv')
SPM_canonicalHRF_P2 = pd.read_csv(SPM_canonicalHRF_P2, sep=';', header=None)
SPM_canonicalHRF_P2 = SPM_canonicalHRF_P2.loc[:, 0]

# load in the SPM canonical HRF - with default peak at 1s (instead of 6s) 
SPM_canonicalHRF_P1 = os.path.join(group_HRF_path, 'SPM_canonicalHRF_shift1.csv')
SPM_canonicalHRF_P1 = pd.read_csv(SPM_canonicalHRF_P1, sep=';', header=None)
SPM_canonicalHRF_P1 = SPM_canonicalHRF_P1.loc[:, 0]

HRF_canonical = SPM_canonicalHRF.loc[:, 0] # seperate them into temporal and dispersion derivs
HRF_temporal = SPM_canonicalHRF.loc[:, 1]
HRF_dispersion = SPM_canonicalHRF.loc[:, 2]

#------------------------------------------------------------------------------
# load in the ROI hrf file -> to add later

ROI_HRF_file=os.path.join(ROI_HRF_path, '1_canonical', 'rsHRF_both_days.csv')
ROI_HRF_file = pd.read_csv(ROI_HRF_file, sep=',')


#------------------------------------------------------------------------------
# load in the BOLD LC onsets according to the rs-HRF analysis 
LC_BOLD_onsets_path=os.path.join(Edrive_ROI_HRF_path, '1_canonical', 'event_onsets')


# specify day
day = ['day1', 'day2']

header_prop = ['subj', 'day', 'prop_invalid']
l_prop = []
l_subject = []
l_day = []


save_dat = []
save_subj = []
save_day = []

timecourse_all=[]
BOLD_onsets_all = []
BOLD_onsets_total = []
deconv_LC  = []
epochs_total_BOLD=[]
epochs_total_pup=[]
epoch_raw_signal_total =make_empty_array(0, baseline_dur+epoch_dur, sF)
#==============================================================================
# Loop over subjects
for subj in range(len(os.listdir(raw_datapath))):

    if subj == 17:
        print('stop')
        continue
    if subj == 18:
        print('stop')
        continue
    if subj == 20:
        print('stop')
        continue
    if subj == 23:
        print('stop')
        continue
    if subj == 28:
        print('stop')
        continue
#            if d == 0:
#                print('stop')
#                continue
    if subj==30:
        continue
    l_pup_int_total=[]
    # Loop over days
    for d in range(len(day)):

        # check for > 25% interopleated data 
        filename = glob.glob(os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], "*.asc"))
        prop = get_invalid_samples(filename[0],sF)
        
        # save % of invalid samples to seperate dataframe 
        l_subject.append(os.listdir(raw_datapath)[subj])
        l_day.append(day[d])
        l_prop.append(prop)
        # combine list
        l_prop_invalid = [l_subject, l_day, l_prop]
        l_prop_invalid=list(map(list, zip(*l_prop_invalid))) #flip   
        l_prop_invalid.insert(0,header_prop)
        l_prop_invalid=pd.DataFrame(l_prop_invalid)
        # Save group .csv file
        #l_prop_invalid.to_csv(os.path.join(statspath,'prop_invalid_samples.csv'))
        
        # Name current file 
        filename = glob.glob(os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], 'PupCor_output',  "MRI*.txt"))
#        if len(filename) < 2:
#             continue
        
    #--------------------------------------------------------------------------
    # Load in smoothed pup data

            
        l_pup_int = []
        
        if len(filename) == 2:
        
            with open(filename[1], 'r') as filehandle:
                for line in filehandle:
                    # remove linebreak which is the last character of the string
                    pup_int = int(float(line[:-1]))
    
                    # add item to the list
                    l_pup_int.append(pup_int)
                    
    #--------------------------------------------------------------------------
    # MAKE  DIRS
            #if processed dir not there, make it
            save_dir_processed_smoothed = os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], 'processed', 'smoothed')
            if not os.path.exists(save_dir_processed_smoothed):
                os.makedirs(save_dir_processed_smoothed)
    
    
#--------------------------------------------------------------------------
# Correct the length to the start and end time of RS

            #take the first 5min (15000 samples) since there is a black screen at the end .
            endtime=nScans*TR*sF
            l_pup_int = l_pup_int[0:endtime]
            
 #--------------------------------------------------------------------------          
# Log the subject number and pupil variance (standard deviation) - will save this later
            save_dat.append(np.std(l_pup_int))
            save_day.append(d)
            save_subj.append(os.listdir(raw_datapath)[subj])
            

 #--------------------------------------------------------------------------
# Apply Band-pass filter?
            


                             
#--------------------------------------------------------------------------
# Calculate the 1st order derivative  of pupil size
            deriv_pup_int = list(np.diff(l_pup_int))
            deriv_pup_int.insert(0, 0)
            deriv_pup_int=np.array(deriv_pup_int)
               
 
# Save data 'orig':
            # save the processed data (1_orig)
            save_smooth_dir = os.path.join(save_dir_processed_smoothed, "1_orig")
            if not os.path.exists(save_smooth_dir):
                os.makedirs(save_smooth_dir)
                
            print('saving file:', filename[1])
        
            np.savetxt(os.path.join(save_smooth_dir, "pupil_dilation.txt"), l_pup_int, delimiter=',')
            np.savetxt(os.path.join(save_smooth_dir, "pupil_derivative.txt"), deriv_pup_int, delimiter=',')


# concatonate the pup timecourse across both days (per subj) 
            l_pup_int_total.append(l_pup_int)

#--------------------------------------------------------------------------
# Shift pupil timecourse back in time by 1s
            #shift backwards by 1 second to account for lag
            l_pup_int_m1s = l_pup_int[sF:15000]  # sF = 50 = 1sec
            zeros = [0] * sF
            l_pup_int_m1s = np.insert(l_pup_int_m1s, 14950, zeros)
    
            deriv_pup_int_m1s = deriv_pup_int[sF:15000]
            deriv_pup_int_m1s = np.insert(deriv_pup_int_m1s, 14950, zeros)

#--------------------------------------------------------------------------
# Do downsampling to TR
       
            # Downsample to TR (2s) 
            pup_int_TR_bins = downsample2bins(l_pup_int_m1s, TR*sF, nScans, 1)
            deriv_pup_int_TR = downsample2bins(deriv_pup_int_m1s, TR*sF, nScans, 1)
            
            l_pup_int_DS_no1sec = downsample2bins(l_pup_int, TR*sF, nScans, 1)
            deriv_pup_int_TR_no1sec = downsample2bins(deriv_pup_int, TR*sF, nScans, 1)
            
            
#--------------------------------------------------------------------------              
   # Center the pupil vectors              

             #center the PUP int TR bin 
            pup_int_TR_bins_center = pup_int_TR_bins-np.average(pup_int_TR_bins)
            
            l_pup_int_DS_no1sec_center = l_pup_int_DS_no1sec-np.average(l_pup_int_DS_no1sec)
            deriv_pup_int_TR_no1sec_center = deriv_pup_int_TR_no1sec-np.average(deriv_pup_int_TR_no1sec)
            
            
            #center the PUP deriv TR bin 
            deriv_pup_int_TR_center = deriv_pup_int_TR-np.average(deriv_pup_int_TR)
            
            
            
            
            
#--------------------------------------------------------------------------
# Shift pupil timecourse back in time by 2s   

            l_pup_int_plus2s = np.insert(l_pup_int_DS_no1sec_center, 0, 0)
            l_pup_int_plus2s = l_pup_int_plus2s[0:150]  # sF = 50 = 1sec
 
#--------------------------------------------------------------------------
# RAW  PUPILS MODELS :           
#--------------------------------------------------------------------------
# Make 5 raw data modes: -2bin, -1bin, 0bin, +1bin, +2bin
            
            # plus 1 TR 
            # size
            pup_int_TR_bins_p1bin=np.insert(pup_int_TR_bins_center, 0,0)
            # deriv
            deriv_pup_int_TR_p1bin=np.insert(deriv_pup_int_TR_center,0,0)
            
            
            # plus 2 TR
             # size
            pup_int_TR_bins_p2bin=np.insert(pup_int_TR_bins_p1bin,0,0)
            # deriv
            deriv_pup_int_TR_p2bin=np.insert(deriv_pup_int_TR_p1bin,0,0)
            
            
            # minus 1 TR
            # size
            pup_int_TR_bins_m1bin=pup_int_TR_bins_center[1:150]
            pup_int_TR_bins_m1bin=np.append(pup_int_TR_bins_m1bin, 0)
            # deriv
            deriv_pup_int_TR_m1bin=deriv_pup_int_TR_center[1:150]
            deriv_pup_int_TR_m1bin=np.append(deriv_pup_int_TR_m1bin, 0)
            
            # minus 2 TR
            # size
            pup_int_TR_bins_m2bin=pup_int_TR_bins_center[2:150]
            pup_int_TR_bins_m2bin=np.append(pup_int_TR_bins_m2bin, 0)
            pup_int_TR_bins_m2bin=np.append(pup_int_TR_bins_m2bin, 0)
            # deriv
            deriv_pup_int_TR_m2bin=deriv_pup_int_TR_center[2:150]
            deriv_pup_int_TR_m2bin=np.append(deriv_pup_int_TR_m2bin, 0)
            deriv_pup_int_TR_m2bin=np.append(deriv_pup_int_TR_m2bin, 0)
            
            
# MAKE DIRECTORIES        
            save_dir_2_smoothed = os.path.join(save_dir_processed_smoothed, '2_downsampled')
            if not os.path.exists(save_dir_2_smoothed):
                os.makedirs(save_dir_2_smoothed)
                
            np.savetxt(os.path.join(save_dir_2_smoothed,"pupil_dilation.txt"), l_pup_int_DS_no1sec_center, delimiter=',')
            np.savetxt(os.path.join(save_dir_2_smoothed,"pupil_derivative.txt"), deriv_pup_int_TR_no1sec_center, delimiter=',')
   
            save_dir_4_smoothed = os.path.join(save_dir_processed_smoothed, '4_raw_pupil_mods')
            if not os.path.exists(save_dir_4_smoothed):
                os.makedirs(save_dir_4_smoothed)
            

            np.savetxt(os.path.join(save_dir_4_smoothed,"pupil_dilation_0bin.txt"), pup_int_TR_bins_center, delimiter=',')
            np.savetxt(os.path.join(save_dir_4_smoothed,"pupil_derivative_0bin.txt"), deriv_pup_int_TR_center, delimiter=',')
            
            
            np.savetxt(os.path.join(save_dir_4_smoothed,"pupil_dilation_p1_bin.txt"), pup_int_TR_bins_p1bin[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_4_smoothed,"pupil_dilation_p2_bin.txt"), pup_int_TR_bins_p2bin[0:150], delimiter=',')

            
            np.savetxt(os.path.join(save_dir_4_smoothed,"pupil_derivative_p1_bin.txt"), deriv_pup_int_TR_p1bin[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_4_smoothed,"pupil_derivative_p2_bin.txt"), deriv_pup_int_TR_p2bin[0:150], delimiter=',')
            
            
            np.savetxt(os.path.join(save_dir_4_smoothed,"pupil_dilation_m1_bin.txt"), pup_int_TR_bins_m1bin[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_4_smoothed,"pupil_dilation_m2_bin.txt"), pup_int_TR_bins_m2bin[0:150], delimiter=',')

            
            np.savetxt(os.path.join(save_dir_4_smoothed,"pupil_derivative_m1_bin.txt"), deriv_pup_int_TR_m1bin[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_4_smoothed,"pupil_derivative_m2_bin.txt"), deriv_pup_int_TR_m2bin[0:150], delimiter=',')
            
 
    
            
            # Apply band-pass filtering to pupil size 
            #BP_pup_size = butter_bandpass_filter(pup_int_TR_bins_center, 0.01, 0.1, 0.5, order=1)



#--------------------------------------------------------------------------                   
#--------------------------------------------------------------------------
# Do convolution
# Load in Indiv HRF: occipital lobe              
            #get occipital ROI HRF from file for this subject
#                HRF_fin = HRF_fin.loc[HRF_fin['subjname'] == os.listdir(datapath)[subj]]
#                del HRF_fin['subjname']
#                HRF_fin = HRF_fin.values.tolist()
#                HRF_subj = HRF_fin[0]
#                HRF_subj=[val - HRF_subj[0] for val in HRF_subj]
#                
#     
#                #do the convolution [Mode ‘same’ returns output of length max(M, N). Boundary effects are still visible.]
#                pup_int_TR_bins_conv = np.convolve(pup_int_TR_bins_center, HRF_subj, mode='full')
#                plt.plot(pup_int_TR_bins_conv)
#                plt.plot(pup_int_TR_bins_center)
#                scaleHRF = np.interp(HRF_subj, (min(HRF_subj), max(HRF_subj)), (min(HRF_SPM32), max(HRF_SPM32)))
#                
#            	#pup_int_TR_bins_conv_scaled = np.convolve(pup_int_TR_bins_center, scaleHRF, mode='full')
#                #pup_int_TR_bins_conv_scaled_MIN2 = np.convolve(pup_int_TR_bins_center, scaleHRF[1:], mode='full')
#                
#                scaleHRF = np.interp(HRF_subj, (min(HRF_subj), max(HRF_subj)), (min(HRF_SPM32), max(HRF_SPM32)))
#
#                
#     # SAVE CONVOLVED W/INDIV HRF DATA                    
#                #make directory: smoothed and unsmoothed dir for convolved data
#                save_dir_3_unsmoothed = os.path.join(l_process_dirs[0], '3_convolved', 'HRF_occ')
#                if not os.path.exists(save_dir_3_unsmoothed):
#                    os.makedirs(save_dir_3_unsmoothed)
#
#                save_dir_3_smoothed = os.path.join(l_process_dirs[1], '3_convolved', 'HRF_occ')
#                if not os.path.exists(save_dir_3_smoothed):
#                    os.makedirs(save_dir_3_smoothed)
#                    


#                    np.savetxt(os.path.join(save_dir_3_smoothed,"pupil_dilation.txt"), pup_int_TR_bins_zscore_conv, delimiter=',')
#                    np.savetxt(os.path.join(save_dir_3_smoothed,"pupil_derivative.txt"), deriv_pup_int_TR_bins_zscore_conv, delimiter=',')
#                
 

 #--------------------------------------------------------------------------
# Do convolution
# Canonical HRF

           #make directory: smoothed and unsmoothed dir for convolved data

            save_dir_3_smoothed_1 = os.path.join(save_dir_processed_smoothed, '3_convolved', 'HRF_canonical')
            if not os.path.exists(save_dir_3_smoothed_1):
                os.makedirs(save_dir_3_smoothed_1)
           
            # temporal
            save_dir_3_smoothed_2 = os.path.join(save_dir_processed_smoothed, '3_convolved', 'HRF_temporal')
            if not os.path.exists(save_dir_3_smoothed_2):
                os.makedirs(save_dir_3_smoothed_2)
          
            # derivative
            save_dir_3_smoothed_3 = os.path.join(save_dir_processed_smoothed, '3_convolved', 'HRF_dispersion')
            if not os.path.exists(save_dir_3_smoothed_3):
                os.makedirs(save_dir_3_smoothed_3)
                
 #Make the canonical HRF       
            # prepare the HRF file to array
            HRF_canonical = np.array(HRF_canonical)
            HRF_temporal = np.array(HRF_temporal)
            HRF_dispersion = np.array(HRF_dispersion)
            
#Make the canonical HR - with peak at 5s
            SPM_canonicalHRF_P5 = np.array(SPM_canonicalHRF_P5)
            
            
 #Convolve with pupil int           
            # Murphy: pup not pushed back 1s  
            # convolve with orig HRF
            pup_int_TR_bins_conv_canonical_murphy = np.convolve(l_pup_int_DS_no1sec_center, HRF_canonical, mode='full')
              # convolve with temp deriv HRF
            pup_int_TR_bins_conv_temporal_murphy = np.convolve(l_pup_int_DS_no1sec_center, HRF_temporal, mode='full')
            # convolve with disper HRF
            pup_int_TR_bins_conv_dispersion_murphy = np.convolve(l_pup_int_DS_no1sec_center, HRF_dispersion, mode='full')

            # TTP: 
            # convolve with orig HRF with peak= 5  (pup pushed back 1s)
            pup_int_TR_bins_conv_canonical = np.convolve(pup_int_TR_bins_center, HRF_canonical, mode='full')
            pup_int_TR_bins_conv_canonical_P5 = np.convolve(pup_int_TR_bins_center, SPM_canonicalHRF_P5, mode='full')
            pup_int_TR_bins_conv_canonical_P4 = np.convolve(pup_int_TR_bins_center, SPM_canonicalHRF_P4, mode='full')
            pup_int_TR_bins_conv_canonical_P3 = np.convolve(pup_int_TR_bins_center, SPM_canonicalHRF_P3, mode='full')
            pup_int_TR_bins_conv_canonical_P2 = np.convolve(pup_int_TR_bins_center, SPM_canonicalHRF_P2, mode='full')
            pup_int_TR_bins_conv_canonical_P1 = np.convolve(pup_int_TR_bins_center, SPM_canonicalHRF_P1, mode='full')
            

        

#Convolve with pupil deriv
            # Murphy: pup not pushed back 1s  
            # convolve with orig HRF
            deriv_pup_int_TR_bins_conv_canonical_murphy = np.convolve(deriv_pup_int_TR_no1sec_center, HRF_canonical, mode='full')
            # convolve with temp deriv HRF
            deriv_pup_int_TR_bins_conv_temporal_murphy = np.convolve(deriv_pup_int_TR_no1sec_center, HRF_temporal, mode='full')
            # convolve with disper HRF
            deriv_pup_int_TR_bins_conv_dispersion_murphy = np.convolve(deriv_pup_int_TR_no1sec_center, HRF_dispersion, mode='full')
            
            # convolve with orig HRF with peak= 5
            deriv_pup_int_TR_bins_conv_canonical = np.convolve(deriv_pup_int_TR_center, HRF_canonical, mode='full')
            deriv_pup_int_TR_bins_conv_canonical_P5 = np.convolve(deriv_pup_int_TR_center, SPM_canonicalHRF_P5, mode='full')
            deriv_pup_int_TR_bins_conv_canonical_P4 = np.convolve(deriv_pup_int_TR_center, SPM_canonicalHRF_P4, mode='full')
            deriv_pup_int_TR_bins_conv_canonical_P3 = np.convolve(deriv_pup_int_TR_center, SPM_canonicalHRF_P3, mode='full')
            deriv_pup_int_TR_bins_conv_canonical_P2 = np.convolve(deriv_pup_int_TR_center, SPM_canonicalHRF_P2, mode='full')
            deriv_pup_int_TR_bins_conv_canonical_P1 = np.convolve(deriv_pup_int_TR_center, SPM_canonicalHRF_P1, mode='full')
            
          
            
            # SAVE FILES   
            #canonical
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_dilation_canonical_murphy.txt"), pup_int_TR_bins_conv_canonical_murphy[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_dilation_canonical.txt"), pup_int_TR_bins_conv_canonical[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_dilation_canonical_P5.txt"), pup_int_TR_bins_conv_canonical_P5[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_dilation_canonical_P4.txt"), pup_int_TR_bins_conv_canonical_P4[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_dilation_canonical_P3.txt"), pup_int_TR_bins_conv_canonical_P3[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_dilation_canonical_P2.txt"), pup_int_TR_bins_conv_canonical_P2[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_dilation_canonical_P1.txt"), pup_int_TR_bins_conv_canonical_P1[0:150], delimiter=',')
           # np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_dilation_canonical_BPF.txt"), pup_int_TR_bins_conv_canonical_BPF[0:150], delimiter=',')
            
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_derivative_canonical_murphy.txt"), deriv_pup_int_TR_bins_conv_canonical_murphy[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_derivative_canonical.txt"), deriv_pup_int_TR_bins_conv_canonical[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_derivative_canonical_P5.txt"), deriv_pup_int_TR_bins_conv_canonical_P5[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_derivative_canonical_P4.txt"), deriv_pup_int_TR_bins_conv_canonical_P4[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_derivative_canonical_P3.txt"), deriv_pup_int_TR_bins_conv_canonical_P3[0:150], delimiter=',') 
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_derivative_canonical_P2.txt"), deriv_pup_int_TR_bins_conv_canonical_P2[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_derivative_canonical_P1.txt"), deriv_pup_int_TR_bins_conv_canonical_P1[0:150], delimiter=',')
            
            
            #save temporal pupil vec
            np.savetxt(os.path.join(save_dir_3_smoothed_2,"pupil_dilation_temporal_murphy.txt"), pup_int_TR_bins_conv_temporal_murphy[0:150], delimiter=',')
           # np.savetxt(os.path.join(save_dir_3_smoothed_2,"pupil_dilation_temporal_BPF.txt"), pup_int_TR_bins_conv_temporal_BPF[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_2,"pupil_derivative_temporal_murphy.txt"), deriv_pup_int_TR_bins_conv_temporal_murphy[0:150], delimiter=',')
                
            #save dispersion pupil vec
            np.savetxt(os.path.join(save_dir_3_smoothed_3,"pupil_dilation_dispersion_murphy.txt"), pup_int_TR_bins_conv_dispersion_murphy[0:150], delimiter=',')
            #np.savetxt(os.path.join(save_dir_3_smoothed_3,"pupil_dilation_dispersion_BPF.txt"), pup_int_TR_bins_conv_dispersion_BPF[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_3,"pupil_derivative_dispersion_murphy.txt"), deriv_pup_int_TR_bins_conv_dispersion_murphy[0:150], delimiter=',')  
            
                
  #--------------------------------------------------------------------------
# Do convolution
# LC HRF    
            #make directory: smoothed and unsmoothed dir for LC convolved data
            save_dir_3_smoothed_1 = os.path.join(save_dir_processed_smoothed, '3_convolved', 'HRF_LC')
            if not os.path.exists(save_dir_3_smoothed_1):
                os.makedirs(save_dir_3_smoothed_1)
          
            # prepare the HRF file to array
            LC_HRF = np.array(ave_LC_HRF)
            LC_HRF=np.ravel(ave_LC_HRF) # make 1dimentional array
            
            #for visuals --> upsample
            f = signal.resample(LC_HRF, 50)
            c = signal.resample(HRF_canonical, 50)
            
#            plt.plot(f, 'y')
#            plt.plot(c, '#FF8C00')
#            
#Convolve with pupil int
            pup_int_TR_bins_conv_LC_HRF = np.convolve(pup_int_TR_bins_center, LC_HRF, mode='full') 
            
            #pup_int_TR_bins_conv_LC_HRF_BPF = np.convolve(BP_pup_size, LC_HRF, mode='full')
            
#Convolve with pupil deriv
            deriv_pup_int_TR_bins_conv_LC_HRF = np.convolve(deriv_pup_int_TR_center, LC_HRF, mode='full')    
            
            
            
              # SAVE FILES   
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_dilation_LC.txt"), pup_int_TR_bins_conv_LC_HRF[0:150], delimiter=',')
           # np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_dilation_LC_BPF.txt"), pup_int_TR_bins_conv_LC_HRF_BPF[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_3_smoothed_1,"pupil_derivative_LC.txt"), deriv_pup_int_TR_bins_conv_LC_HRF[0:150], delimiter=',')                                      
            
            
            
#--------------------------------------------------------------------------
# rsHRF ROI CONVOLUTION: VTA, MR, DR, LC, SN
#--------------------------------------------------------------------------
            
            #define the HRF file    
            roi_rsHRF = ROI_HRF_file
            
            subj_code = os.listdir(raw_datapath)[subj] # get subject code 
            
            # get averabe of MR, LC, DR roi HRFs (replace odd HRFs)
            subj_HRF = roi_rsHRF[roi_rsHRF['subj'] == subj_code]   #get persons HRF for each ROI 
            
            DR_HRF = np.array(list(range(1,len(subj_HRF)+1)))
            LC_HRF = np.array(list(range(1,len(subj_HRF)+1)))
            MR_HRF = np.array(list(range(1,len(subj_HRF)+1)))
            
            for s in set(roi_rsHRF['subj']):
                subj_hrf=roi_rsHRF[roi_rsHRF['subj'] == s]
                DR_HRF=np.vstack([DR_HRF,subj_hrf['DR_roi']])
                LC_HRF=np.vstack([LC_HRF,subj_hrf['LC_roi']])
                MR_HRF=np.vstack([MR_HRF,subj_hrf['MR_roi']])
                
            #calc mean HRF for these ROIs
            DR_HRF_mean = np.nanmean(DR_HRF[1:],axis=0)
            LC_HRF_mean = np.nanmean(LC_HRF[1:],axis=0)
            MR_HRF_mean = np.nanmean(MR_HRF[1:],axis=0)
            
            
            # convolve the centered time adjusted pupil timecourse with each HRF from ROI 
            #subjs: 006, 023, 064
            if subj == 4:
                MR_hrf=signal.resample(MR_HRF_mean[:200], 16)
            else:   
                MR_hrf=signal.resample(subj_HRF['MR_roi'][:200], 16)
                
            if subj == 19:
                LC_hrf=signal.resample(LC_HRF_mean[:200], 16)
            else:
                LC_hrf=signal.resample(subj_HRF['LC_roi'][:200], 16)
            if subj == 51:
                DR_hrf=signal.resample(DR_HRF_mean[:200], 16)
            else:
                DR_hrf=signal.resample(subj_HRF['DR_roi'][:200], 16)
        
            
            # pup size
             # downsample HRF (so convolution works)
            pup_int_conv_DR_roi = np.convolve(pup_int_TR_bins_center,DR_hrf, mode='full')
            
            pup_int_conv_MR_roi = np.convolve(pup_int_TR_bins_center, MR_hrf, mode='full')
            
            pup_int_conv_LC_roi = np.convolve(pup_int_TR_bins_center, LC_hrf, mode='full')
            
            VTA_hrf=signal.resample(subj_HRF['VTA_roi'][:200], 16)
            pup_int_conv_VTA_roi = np.convolve(pup_int_TR_bins_center, VTA_hrf, mode='full')
            
            SN_hrf=signal.resample(subj_HRF['SN_roi'][:200], 16)
            pup_int_conv_SN_roi = np.convolve(pup_int_TR_bins_center, SN_hrf, mode='full')
            
            ACC_hrf=signal.resample(subj_HRF['ACC_roi'][:200], 16)
            pup_int_conv_ACC_roi = np.convolve(pup_int_TR_bins_center, ACC_hrf, mode='full')
            
            OCC_hrf=signal.resample(subj_HRF['OCC_roi'][:200], 16)
            pup_int_conv_OCC_roi = np.convolve(pup_int_TR_bins_center, OCC_hrf, mode='full')
            
            BF_sept_hrf=signal.resample(subj_HRF['BF_sept_roi'][:200], 16)
            pup_int_conv_BF_sept_roi = np.convolve(pup_int_TR_bins_center, BF_sept_hrf, mode='full')
            
            BF_subl_hrf=signal.resample(subj_HRF['BF_subl_roi'][:200], 16)
            pup_int_conv_BF_subl_roi = np.convolve(pup_int_TR_bins_center, BF_subl_hrf, mode='full')
            
             # pup deriv
            deriv_pup_int_conv_DR_roi = np.convolve(deriv_pup_int_TR_center, DR_hrf, mode='full')
            
            deriv_pup_int_conv_MR_roi = np.convolve(deriv_pup_int_TR_center, MR_hrf, mode='full')
            
            deriv_pup_int_conv_LC_roi = np.convolve(deriv_pup_int_TR_center, LC_hrf, mode='full')
            
            deriv_pup_int_conv_VTA_roi = np.convolve(deriv_pup_int_TR_center, VTA_hrf, mode='full')
            
            deriv_pup_int_conv_SN_roi = np.convolve(deriv_pup_int_TR_center, SN_hrf, mode='full')
            
            deriv_pup_int_conv_ACC_roi = np.convolve(deriv_pup_int_TR_center, ACC_hrf, mode='full')

            deriv_pup_int_conv_OCC_roi = np.convolve(deriv_pup_int_TR_center, OCC_hrf, mode='full')
            
            deriv_pup_int_conv_BF_sept_roi = np.convolve(deriv_pup_int_TR_center, BF_sept_hrf, mode='full')
            
            deriv_pup_int_conv_BF_subl_roi = np.convolve(deriv_pup_int_TR_center, BF_subl_hrf, mode='full')
            
            
            save_dir_5_smoothed_1 = os.path.join(save_dir_processed_smoothed, '3_convolved', 'HRF_ROIs')
            if not os.path.exists(save_dir_5_smoothed_1):
                os.makedirs(save_dir_5_smoothed_1)
                
                
            # SAVE FILES   
            # pup size
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_dilation_DR_roi.txt"), pup_int_conv_DR_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_dilation_MR_roi.txt"), pup_int_conv_MR_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_dilation_LC_roi.txt"), pup_int_conv_LC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_dilation_VTA_roi.txt"), pup_int_conv_VTA_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_dilation_SN_roi.txt"), pup_int_conv_SN_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_dilation_ACC_roi.txt"), pup_int_conv_ACC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_dilation_OCC_roi.txt"), pup_int_conv_OCC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_dilation_BF_sept_roi.txt"), pup_int_conv_BF_sept_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_dilation_BF_subl_roi.txt"), pup_int_conv_BF_subl_roi[0:150], delimiter=',')
            
            # pup deriv
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_derivative_DR_roi.txt"), deriv_pup_int_conv_DR_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_derivative_MR_roi.txt"), deriv_pup_int_conv_MR_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_derivative_LC_roi.txt"), deriv_pup_int_conv_LC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_derivative_VTA_roi.txt"), deriv_pup_int_conv_VTA_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_derivative_SN_roi.txt"), deriv_pup_int_conv_SN_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_derivative_ACC_roi.txt"), deriv_pup_int_conv_ACC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_derivative_OCC_roi.txt"), deriv_pup_int_conv_OCC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_derivative_BF_sept_roi.txt"), deriv_pup_int_conv_BF_sept_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(save_dir_5_smoothed_1,"pupil_derivative_BF_subl_roi.txt"), deriv_pup_int_conv_BF_subl_roi[0:150], delimiter=',')
            
                                        
    unpack_l_pup_int_total = [j for i in l_pup_int_total for j in i]     
     #------------------------------------------------------------------------------
    # run point-process events on the pupil timecourse
    threshold=1
    #centre pupil timecourse 
    centered_concat_timecourse = unpack_l_pup_int_total-np.average(unpack_l_pup_int_total)
    centered_concat_timecourse = np.array(centered_concat_timecourse)
    
    # concatonate pupil timecourse per subject
    centered_concat_timecourse_total = {'subject': os.listdir(raw_datapath)[subj],
                                        't': list(range(0,3001)),
                                        'pupil_timecourse': centered_concat_timecourse}
    timecourse_all.append(centered_concat_timecourse_total)
    # run point-process events on the pupil timecourse
    onsets_pup = get_point_processEvents(unpack_l_pup_int_total, 2, sF)


    rf = nd.ResponseFitter(input_signal=centered_concat_timecourse,
                    sample_rate=sF)
    
    # point-process based on events in pup
    cue_epochs = rf.get_epochs(onsets=onsets_pup,
                   interval=[0, 10.5])
    
#    
#    
#    ## plots pupil response ! =================================================
##    palette = sns.color_palette('Set1')
##    
#    epochs_total_pup.append(cue_epochs.mean())
#    epochs = np.array(epochs_total_pup)

#    
#    # plot averave epochs 
#    epochs_mean = np.nanmean(epochs[:],axis=0)
#    epochs_SEM = stats.sem(epochs[:], axis=0, nan_policy= 'omit')
#    epochs_SEM1 = epochs_mean - epochs_SEM 
#    epochs_SEM2 = epochs_mean + epochs_SEM
#    x = np.array(list(range(0,int((0--10)*sF)+1)))
#    
#    # plot: 
#    plt.plot(epochs_mean) 
#    plt.fill_between(x, epochs_SEM1,epochs_SEM2, alpha = 0.4)
#    plt.axvline(x=sF, color='black', linestyle='--')
#    plt.xlabel("Time [sec]")
#    plt.ylabel("Pupil diameter [a.u]")
#    plt.title('Epoched average of psuedo-event-locked response')
#    
#    tick_locs = np.arange(0, 501, 100).tolist()
#    tick_lbls = np.arange(0, 12, 2).tolist()
#    plt.xticks(tick_locs,tick_lbls)
#    plt.figtext(0.1, 0.05, '----- pupil crosses 2SD threshold', fontsize=8)
#
#
#    
#    cue_epochs.T.plot(c=palette[0], alpha=.5, ls='--', legend=False)
#    cue_epochs.mean().plot(c=palette[0], lw=2, alpha=1.0)
#    plt.axvline(x=sF, color='black', linestyle='--')
#    sns.despine()
#    plt.xlabel('Time (s)')
#    plt.title('Epoched average of cue-locked response')
#    plt.gcf().set_size_inches(8, 3)
#    

    
    # ROI definion: 
    rois = ['LC', 'VTA', 'SN', 'DR', 'MR', 'BF_subl', 'ACC', 'OCC']
    
    for c_roi in rois:

        #------------------------------------------------------------------------------
        # Run similar point-process on the pupil, but this time using the LC BOLD events (rs-HRF analysis)

        BOLD_onsets = os.path.join(LC_BOLD_onsets_path, (os.listdir(raw_datapath)[subj] + '_' + c_roi + '_rsHRF_event_onsets_both_days.csv'))
        BOLD_onsets = pd.read_csv(BOLD_onsets, sep=',')
        BOLD_onsets = BOLD_onsets[c_roi]
        BOLD_onsets = [x*2 for x in BOLD_onsets]   # this is because the BOLD onsets are measure in TR - need to conver to seconds
        
    
    
    
        #pd.DataFrame(data=centered_concat_timecourse_total)
        # concatonate BOLD onsets per subject
        BOLD_onsets_total = {'subject': os.listdir(raw_datapath)[subj],
                             'trial_type': c_roi, 
                             'BOLD_onset': BOLD_onsets}
        
        
        BOLD_onsets_all.append(BOLD_onsets_total)
BOLD_onsets_all_dat = pd.DataFrame(data=BOLD_onsets_all)
timecourse_all_dat = pd.DataFrame(data=timecourse_all)

filepath = 'pupil_timecourse.csv' 
timecourse_all_dat.to_csv(filepath)
g_model = GroupResponseFitter(timecourse_all_dat,
                              BOLD_onsets_all_dat,
                              input_sample_rate=50,
                              concatenate_runs=False)
g_model.add_event('LC',
                  basis_set='fourier',
                  n_regressors=9,
                  interval=[0, 8])

        # get average (centered) pupil timcourse for these events (+ 10s )
#        epoch_sig = extract_ave_epoch(centered_concat_timecourse, BOLD_onsets, baseline_dur, epoch_dur, sF)
#        epoch_sig = np.nanmean(epoch_sig[1:],axis=0)  #averege the subject signal 
#        
#        # plots the indiv pupil signal in seperate plots 
#        plt.plot(epoch_sig)
#        plt.savefig('C:\\Users\\lloydb\\surfdrive\\ExperimentData\\NYU_RS_LC\\other\\manuscript\\figures\\fresh_pup_figs\\indiv_participants\\' + c_roi + '\\' 
#                    + os.listdir(raw_datapath)[subj] + '_ave_pupresponse_to_BOLD.png', format='png')
#        plt.close()
#        
#        epoch_raw_signal_total = np.vstack([epoch_raw_signal_total, epoch_sig])  # combine subjects 
#
#        np.array(epochs_total_BOLD)
#        
#        
#        # deconvolve pup based on LC events
#        cue_epochs = rf.get_epochs(onsets=BOLD_onsets,
#                       interval=[0, 11])
        
        
#        len_epoch=6
#        rf.add_event(event_name='BOLD spike',
#             onsets=BOLD_onsets,
#             interval=[0,len_epoch])
#        
#        # save indiv deconv response 
#        rf.fit()
#        x = rf.get_timecourses()
#        rf.plot_timecourses()
#        plt.suptitle('Linear deconvolution using GLM and FIR')
#        plt.title('')
#        plt.savefig('C:\\Users\\lloydb\\surfdrive\\ExperimentData\\NYU_RS_LC\\other\\manuscript\\figures\\fresh_pup_figs\\indiv_deconv\\' + c_roi + '\\' 
#            + os.listdir(raw_datapath)[subj] + '_ave_deconv_pup_BOLD.png', format='png')
#        
#        plt.close()
#        deconv_LC.append(x)
        
#        epochs_total_BOLD.append(cue_epochs.mean())
        
        print('moving along subject-wise')
        
## first plot deconvolved respons! 
#epochs_plot = np.array(epochs_total_BOLD)
#
## plot deconvolved events!
#palette = sns.color_palette('Set1')
#FS = 20
#FS_tick = 15
#
## plot average signal ! 
#epochs_mean = np.nanmean(epochs_plot[:],axis=0)
#epochs_SEM = stats.sem(epochs_plot[:], axis=0, nan_policy= 'omit')
#epochs_SEM1 = epochs_mean - epochs_SEM 
#epochs_SEM2 = epochs_mean + epochs_SEM
#x = np.array(list(range(0,int((0--11)*sF)+1)))
#
## plot:
#plt.rcParams['font.size'] = FS_tick
#plt.plot(epochs_mean, color  = 'black') 
#plt.fill_between(x, epochs_SEM1,epochs_SEM2)
#plt.xlabel("Time [sec]", fontsize = FS)
#plt.ylabel("Pupil diameter [a.u]", fontsize = FS)
#plt.title(c_roi, fontsize = FS)
#tick_locs = np.arange(0, 551, 50).tolist()
#tick_lbls = np.arange(0, 11, 1).tolist()
#plt.xticks(tick_locs,tick_lbls, fontsize=FS_tick)
#
#
#
## get arrays in correct way
#epochs_mean=np.nanmean(deconv_LC[:],axis=0)
#epoch_sem=stats.sem(deconv_LC[:], axis=0, nan_policy= 'omit')
#mean_flattened = [] 
#for sublist in epochs_mean: 
#    for val in sublist: 
#        mean_flattened.append(val) 
#sem_flattened = [] 
#for sublist in epoch_sem: 
#    for val in sublist: 
#        sem_flattened.append(val)         
#        
#
#
## plot ave deconv response! 
#palette = sns.color_palette('Set1')
#FS = 20
#FS_tick = 15
#sF_devonv = 1000
#
## plot averave epochs 
#epochs_mean = mean_flattened
#epochs_SEM = sem_flattened
#epochs_SEM1 = np.array(epochs_mean) - np.array(epochs_SEM) 
#epochs_SEM2 = np.array(epochs_mean) + np.array(epochs_SEM) 
#x = np.array(list(range(1,int((0--len_epoch)*sF_devonv)+1)))
#
## plot:
#plt.rcParams['font.size'] = FS_tick
#plt.plot(epochs_mean, color  = 'black') 
#plt.fill_between(x, epochs_SEM1,epochs_SEM2)
#plt.xlabel("Time [sec]", fontsize = FS)
#plt.ylabel("Pupil diameter [a.u]", fontsize = FS)
#plt.title(c_roi, fontsize = FS)
#tick_locs = np.arange(0, len_epoch*sF_devonv, sF_devonv).tolist()
#tick_lbls = np.arange(0, len_epoch, 1).tolist()
#plt.xticks(tick_locs,tick_lbls, fontsize=FS_tick)
#
#plt.savefig('C:\\Users\\lloydb\\surfdrive\\ExperimentData\\NYU_RS_LC\\other\\manuscript\\figures\\fresh_pup_figs\\average_deconv_BOLD_ROIS\\'+c_roi+'\\' + c_roi + '_ave_deconv_pup.eps', format='eps')

#plt.axvline(x=0, color='black', linestyle='--')
#plt.figtext(0.1, 0.04, '----- VTA BOLD crosses 1SD threshold', fontsize=8)



#        cue_epochs.T.plot(c=palette[0], alpha=.5, ls='--', legend=False)
#        cue_epochs.mean().plot(c=palette[0], lw=2, alpha=1.0)
#        plt.axvline(x=sF, color='black', linestyle='--')
#        sns.despine()
#        plt.xlabel('Time (s)')
#        plt.title('Epoched average of cue-locked response')
#        plt.gcf().set_size_inches(8, 3)
#                
    
               
#--------------------------------------------------------------------------
# Match events between the three onsets: pupil, LC, VTA    
total_onsets = []
condition = []
for i, x in enumerate(onsets_pup):  
    # onsets_df['pupil'][round(x)] = 1
     total_onsets.append(round(x))
     condition.append('pup')
for i, x in enumerate(onsets_LC):  
    # onsets_df['LC'][round(x)] = 1
     total_onsets.append(x)
     condition.append('LC')
for i, x in enumerate(onsets_VTA):  
     #onsets_df['VTA'][round(x)] = 1
     total_onsets.append(x)
     condition.append('VTA')
  
   # total_onsets.insert(0, "Onset")
   # condition.insert(0, 'Condition')
total_lists= {'Onset':total_onsets,'Condition':condition}
df = pd.DataFrame(data=total_lists)
# plot the onset timepoints across conditions 
p = sns.catplot(x="Onset", y="Condition", data=df)
p.set_xlabel("X-Axis", fontsize = 20)
p.set_ylabel("Y-Axis", fontsize = 20)
            
            
            
    
save_all = [save_subj, save_day, save_dat]
save_all=list(map(list, zip(*save_all))) #flip    

save_all_header = ['subject', 'day', 'SD pupil']
    #save
#save_all.insert(0,save_all_header)
#savetxt(save_all, "D:\\NYU_RS_LC\\stats\\pupil\\standard_dev_pupil.txt")




#fig, axs = plt.plot()
#plt.rc('font', size=22)   
##fig.suptitle('pupil size')
#plt.plot.rc('axes', labelsize=16)    # fontsize of the x and y labels
#plt.plot.rc('xtick', labelsize=16)    # fontsize of the tick labels
#plt.plot.rc('ytick', labelsize=16)
plt.figure(figsize=(11, 3))
plt.rc('font', size=32) 
#plt.rc('legend',fontsize=32)  
plt.plot(l_pup_int_DS_no1sec_center[0:149],'tab:blue', linewidth=1) # downsampled pup, centered, no 1-sec shift \
plt.plot(pup_int_TR_bins_conv_canonical_P1[0:149],'tab:red', linewidth=1)   # downsampled pup, +2s, centered 
plt.plot(pup_int_TR_bins_conv_canonical[0:149],'tab:grey', linewidth=1)   # convolved with 1s HRF
plt.plot(sinewave)
plt.ylabel('pupil size')
#plt.plot(0, 0, "blue", label="no lag")
#plt.plot(0, 0, "red", label = "+2s lag")
#plt.plot(0, 0, "green", label = "conv TTP=1")
#plt.legend()

start_time = 0
end_time = 150
sample_rate = 1
time = np.arange(start_time, end_time, 1/sample_rate)
theta = 0
frequency = 0.1
amplitude = 500
sinewave = amplitude * np.sin(2 * np.pi * frequency * time + theta)
plt.plot(sinewave)


plt.savefig('overlap_plots_new.eps', format='eps')

# plots: 
fig, axs = plt.subplots(3)
plt.rc('font', size=22)   
#fig.suptitle('pupil size')
plt.rc('axes', labelsize=22)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=22)    # fontsize of the tick labels
plt.rc('ytick', labelsize=22)
axs[0].plot(l_pup_int,'tab:blue') # raw pup
axs[1].plot(pup_int_TR_bins_center,'tab:red')   # downsampled
axs[2].plot(pup_int_conv_LC_roi[0:150],'tab:green')   # convolved 

plt.sca(axs[0])
plt.xticks(np.arange(0, 16000, step=5000))
plt.yticks([0, 1000, 2000, 3000])

plt.sca(axs[1])
plt.xticks(np.arange(0, 160, step=50))
plt.yticks([-1000, 0, 2000])

plt.sca(axs[2])
plt.xticks(np.arange(0, 160, step=50))
plt.yticks([-1000, 0, 2000])
plt.savefig('process_pupsize_plots.eps', format='eps')

# plots: 
fig, axs = plt.subplots(3)
plt.rc('font', size=22) 
#fig.suptitle('pupil derivative')
plt.rc('axes', labelsize=22)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=22)    # fontsize of the tick labels
plt.rc('ytick', labelsize=22)
axs[0].plot(deriv_pup_int,'tab:blue')  # pup deriv
axs[1].plot(deriv_pup_int_TR_center,'tab:red')# downsampled
axs[2].plot(deriv_pup_int_conv_LC_roi[0:150],'tab:green') # downsampled



#zscore_pup = stats.zscore(pup_int_TR_bins)
#zscore_pup=list(zscore_pup)
#zscore_pup_conv = np.convolve(zscore_pup, SPM_canonicalHRF_P2, mode='full')
#
#lc_data_029 = 'E:\\NYU_RS_LC\\stats\\rsHRF\\MRI_FCWML029\\concat\\1_canonical\\LC_rawdata_zscore.txt'
#lc_data_029 = pd.read_csv(lc_data_029, sep=';', header=None)
#plt.plot(lc_data_029)
#plt.plot(zscore_pup_conv)
#
