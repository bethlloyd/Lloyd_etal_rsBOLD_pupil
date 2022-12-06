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

## ===== SET THIS VARIABLE! ====== ## 
save_data = False
do_rsHRF_data = True

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
group_HRF_path = os.path.join(homepath, 'stats', 'fMRI', '4a_TTP_analysis_estimation')
# to add later:
ROI_HRF_path ="D:\\NYU_RS_LC\\stats\\R_csv_files"
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

ROI_HRF_file=os.path.join(ROI_HRF_path, 'rsHRF_extracted_HRFs_both_days.csv')
ROI_HRF_file = pd.read_csv(ROI_HRF_file, sep=',')

# settings------------------
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

    # Loop over days
    for d in range(len(day)):
        if (os.listdir(raw_datapath)[subj] == 'MRI_FCWML025') and (d == 0):
            continue
        if (os.listdir(raw_datapath)[subj] == 'MRI_FCWML036') and (d == 1):  # no pupil file for this subj
            print(f'skipping subject {os.listdir(raw_datapath)[subj]} day {d+1}')
            continue
        if (os.listdir(raw_datapath)[subj] == 'xMRI_FCWML044'):
            continue
        # check for > 25% interopleated data 
        filename = glob.glob(os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], '1_rawdat', "*.asc"))
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
        l_prop_invalid.to_csv(os.path.join(homepath, 'stats', 'R_csv_files', 'pupil_orig_prop_invalid_samples.csv'))
        
        # Name current file 
        filename = glob.glob(os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], '1_rawdat',  "MRI*smth_int_pup.txt"))
        if len(filename) < 1:
             print(f'No smoothed pup file for subject {os.listdir(raw_datapath)[subj]}')
        
    #--------------------------------------------------------------------------
    # Load in smoothed pup data
        l_pup_int = []
        
        with open(filename[0], 'r') as filehandle:
            for line in filehandle:
                # remove linebreak which is the last character of the string
                pup_int = int(float(line[:-1]))

                # add item to the list
                l_pup_int.append(pup_int)
                
        print(f'processing pupil for subject: {os.listdir(raw_datapath)[subj]}')
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
# 1. Calculate the 1st order derivative  of pupil size
        deriv_pup_int = list(np.diff(l_pup_int))
        deriv_pup_int.insert(0, 0)
        deriv_pup_int=np.array(deriv_pup_int)
         
        
        if save_data:
            # Save data: pupil dilation + pupil derivative --> '2_dat':
            save_2_dat_dir = os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], "2_dat")
            np.savetxt(os.path.join(save_2_dat_dir, "pupil_dilation.txt"), l_pup_int, delimiter=',')
            np.savetxt(os.path.join(save_2_dat_dir, "pupil_derivative.txt"), deriv_pup_int, delimiter=',')


#--------------------------------------------------------------------------
# 2. Downsample pupil timeseries to TR(2s) --> pupil vector used for the cross-corr
        l_pup_int_DS_no1sec = downsample2bins(l_pup_int, TR*sF, nScans, 1)
        deriv_pup_int_TR_no1sec = downsample2bins(deriv_pup_int, TR*sF, nScans, 1)


#--------------------------------------------------------------------------
# 3. Shift pupil timecourse back in time by 1s and downsample --> pupil vector used for HRF analyses 
        l_pup_int_m1s = l_pup_int[sF:15000]  # sF = 50 = 1sec
        zeros = [0] * sF
        l_pup_int_m1s = np.insert(l_pup_int_m1s, 14950, zeros)
        
        deriv_pup_int_m1s = deriv_pup_int[sF:15000]
        deriv_pup_int_m1s = np.insert(deriv_pup_int_m1s, 14950, zeros)

        pup_int_TR_bins = downsample2bins(l_pup_int_m1s, TR*sF, nScans, 1)
        deriv_pup_int_TR = downsample2bins(deriv_pup_int_m1s, TR*sF, nScans, 1)
        
        
# Save data: pupil dilation + pupil derivative --> '2_dat':
        save_3_downsampled_dir = os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], "3_downsampled")
        if save_data:
            np.savetxt(os.path.join(save_3_downsampled_dir,"pupil_dilation_1s_shifted.txt"), pup_int_TR_bins, delimiter=',')
            np.savetxt(os.path.join(save_3_downsampled_dir,"pupil_derivative_1s_shifted.txt"), deriv_pup_int_TR, delimiter=',')
#--------------------------------------------------------------------------              
# 4. Center the pupil vectors              

        l_pup_int_DS_no1sec_center = l_pup_int_DS_no1sec-np.average(l_pup_int_DS_no1sec)
        deriv_pup_int_TR_no1sec_center = deriv_pup_int_TR_no1sec-np.average(deriv_pup_int_TR_no1sec)

        pup_int_TR_bins_center = pup_int_TR_bins-np.average(pup_int_TR_bins)
        deriv_pup_int_TR_center = deriv_pup_int_TR-np.average(deriv_pup_int_TR)



# Save pupil data: downsampled + centered and not shifted by 1 s --> 3_downsampled
        if save_data:
            save_3_downsampled_dir = os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], "3_downsampled")
            np.savetxt(os.path.join(save_3_downsampled_dir,"pupil_dilation.txt"), l_pup_int_DS_no1sec_center, delimiter=',')
            np.savetxt(os.path.join(save_3_downsampled_dir,"pupil_derivative.txt"), deriv_pup_int_TR_no1sec_center, delimiter=',')
        
        
 #--------------------------------------------------------------------------
# 5. Do convolution

# 5a. Murphy replication: canonical + temporal + dispersion: pup not pushed back 1s  
        HRF_canonical = np.array(HRF_canonical)
        HRF_temporal = np.array(HRF_temporal)
        HRF_dispersion = np.array(HRF_dispersion)
        # convolve:
        pup_int_TR_bins_conv_canonical_murphy = np.convolve(l_pup_int_DS_no1sec_center, HRF_canonical, mode='full')
        pup_int_TR_bins_conv_temporal_murphy = np.convolve(l_pup_int_DS_no1sec_center, HRF_temporal, mode='full')
        pup_int_TR_bins_conv_dispersion_murphy = np.convolve(l_pup_int_DS_no1sec_center, HRF_dispersion, mode='full')
        
        # save Murphy replication 
        if save_data:
            save_4_convolved_Murphy_dir = os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], "4_convolved", '1_Murphy_replication')
            if not os.path.exists(save_4_convolved_Murphy_dir):
                os.makedirs(save_4_convolved_Murphy_dir)
            np.savetxt(os.path.join(save_4_convolved_Murphy_dir,"pupil_dilation_canonical.txt"), pup_int_TR_bins_conv_canonical_murphy[0:150], delimiter=',')
            np.savetxt(os.path.join(save_4_convolved_Murphy_dir,"pupil_dilation_temporal.txt"), pup_int_TR_bins_conv_temporal_murphy[0:150], delimiter=',')
            np.savetxt(os.path.join(save_4_convolved_Murphy_dir,"pupil_dilation_dispersion.txt"), pup_int_TR_bins_conv_dispersion_murphy[0:150], delimiter=',')
            
        
# 5b. Schneider replicaton: canonical HRF: pup IS shifted back in time by 1s
        pup_int_TR_bins_conv_canonical_schneider = np.convolve(pup_int_TR_bins_center, HRF_canonical, mode='full')
        # save Schneider replication 
        if save_data:
            save_4_convolved_Schneider_dir = os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], "4_convolved", '2_Schneider_replication')
            if not os.path.exists(save_4_convolved_Schneider_dir):
                os.makedirs(save_4_convolved_Schneider_dir)
            np.savetxt(os.path.join(save_4_convolved_Schneider_dir,"pupil_dilation.txt"), pup_int_TR_bins_conv_canonical_schneider[0:150], delimiter=',')


# 5c. rs-HRF convolution 
        if do_rsHRF_data:
            roi_rsHRF = ROI_HRF_file
            subj_code = os.listdir(raw_datapath)[subj] # get subject code 
            
            # get averabe of MR, LC, DR roi HRFs (replace odd HRFs)
            subj_HRF = roi_rsHRF[roi_rsHRF['subj'] == subj_code]   #get persons HRF for each ROI 

            # convolve ROI HRFs with pupil size:
   
            LC_hrf=signal.resample(subj_HRF['LC'][:200], 16)
            pup_int_conv_LC_roi = np.convolve(pup_int_TR_bins_center, LC_hrf, mode='full')
            
            VTA_hrf=signal.resample(subj_HRF['VTA'][:200], 16)
            pup_int_conv_VTA_roi = np.convolve(pup_int_TR_bins_center, VTA_hrf, mode='full')
            
            SN_hrf=signal.resample(subj_HRF['SN'][:200], 16)
            pup_int_conv_SN_roi = np.convolve(pup_int_TR_bins_center, SN_hrf, mode='full')
            
            DR_hrf=signal.resample(subj_HRF['DR'][:200], 16)
            pup_int_conv_DR_roi = np.convolve(pup_int_TR_bins_center,DR_hrf, mode='full')
            
            MR_hrf=signal.resample(subj_HRF['MR'][:200], 16)
            pup_int_conv_MR_roi = np.convolve(pup_int_TR_bins_center, MR_hrf, mode='full')
            
            BF_hrf=signal.resample(subj_HRF['BF'][:200], 16)
            pup_int_conv_BF_roi = np.convolve(pup_int_TR_bins_center, BF_hrf, mode='full')
            
            ACC_hrf=signal.resample(subj_HRF['ACC'][:200], 16)
            pup_int_conv_ACC_roi = np.convolve(pup_int_TR_bins_center, ACC_hrf, mode='full')
            
            OCC_hrf=signal.resample(subj_HRF['OCC'][:200], 16)
            pup_int_conv_OCC_roi = np.convolve(pup_int_TR_bins_center, OCC_hrf, mode='full')
            
            Pons_hrf=signal.resample(subj_HRF['Pons'][:200], 16)
            pup_int_conv_Pons_roi = np.convolve(pup_int_TR_bins_center, Pons_hrf, mode='full')
            
            
             # pup deriv
            deriv_pup_int_conv_DR_roi = np.convolve(deriv_pup_int_TR_center, DR_hrf, mode='full')
            deriv_pup_int_conv_MR_roi = np.convolve(deriv_pup_int_TR_center, MR_hrf, mode='full')
            deriv_pup_int_conv_LC_roi = np.convolve(deriv_pup_int_TR_center, LC_hrf, mode='full')
            deriv_pup_int_conv_VTA_roi = np.convolve(deriv_pup_int_TR_center, VTA_hrf, mode='full')
            deriv_pup_int_conv_SN_roi = np.convolve(deriv_pup_int_TR_center, SN_hrf, mode='full')
            deriv_pup_int_conv_ACC_roi = np.convolve(deriv_pup_int_TR_center, ACC_hrf, mode='full')
            deriv_pup_int_conv_OCC_roi = np.convolve(deriv_pup_int_TR_center, OCC_hrf, mode='full')   
            deriv_pup_int_conv_BF_roi = np.convolve(deriv_pup_int_TR_center, BF_hrf, mode='full')
            deriv_pup_int_conv_Pons_roi = np.convolve(deriv_pup_int_TR_center, Pons_hrf, mode='full')
            
            
            
            save_4_convolved_rsHRF_dir = os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], "4_convolved", '3_rsHRF_analysis')
            LC_path = os.path.join(save_4_convolved_rsHRF_dir, 'LC')                
            if not os.path.exists(LC_path):
                os.makedirs(LC_path)
            VTA_path = os.path.join(save_4_convolved_rsHRF_dir, 'VTA')                
            if not os.path.exists(VTA_path):
                os.makedirs(VTA_path)
            SN_path = os.path.join(save_4_convolved_rsHRF_dir, 'SN')                
            if not os.path.exists(SN_path):
                os.makedirs(SN_path)
            DR_path = os.path.join(save_4_convolved_rsHRF_dir, 'DR')                
            if not os.path.exists(DR_path):
                os.makedirs(DR_path)
            MR_path = os.path.join(save_4_convolved_rsHRF_dir, 'MR')                
            if not os.path.exists(MR_path):
                os.makedirs(MR_path)
            BF_path = os.path.join(save_4_convolved_rsHRF_dir, 'BF')                
            if not os.path.exists(BF_path):
                os.makedirs(BF_path)
            ACC_path = os.path.join(save_4_convolved_rsHRF_dir, 'ACC')                
            if not os.path.exists(ACC_path):
                os.makedirs(ACC_path)
            OCC_path = os.path.join(save_4_convolved_rsHRF_dir, 'OCC')                
            if not os.path.exists(OCC_path):
                os.makedirs(OCC_path)
            Pons_path = os.path.join(save_4_convolved_rsHRF_dir, 'Pons')                
            if not os.path.exists(Pons_path):
                os.makedirs(Pons_path)
            
                
            # SAVE FILES   
            np.savetxt(os.path.join(DR_path, "pupil_dilation.txt"), pup_int_conv_DR_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(MR_path, "pupil_dilation.txt"), pup_int_conv_MR_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(LC_path, "pupil_dilation.txt"), pup_int_conv_LC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(VTA_path, "pupil_dilation.txt"), pup_int_conv_VTA_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(SN_path, "pupil_dilation.txt"), pup_int_conv_SN_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(ACC_path, "pupil_dilation.txt"), pup_int_conv_ACC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(OCC_path, "pupil_dilation.txt"), pup_int_conv_OCC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(BF_path, "pupil_dilation.txt"), pup_int_conv_BF_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(Pons_path, "pupil_dilation.txt"), pup_int_conv_Pons_roi[0:150], delimiter=',')
            
            # pup deriv
            np.savetxt(os.path.join(DR_path,"pupil_derivative.txt"), deriv_pup_int_conv_DR_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(MR_path,"pupil_derivative.txt"), deriv_pup_int_conv_MR_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(LC_path,"pupil_derivative.txt"), deriv_pup_int_conv_LC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(VTA_path,"pupil_derivative.txt"), deriv_pup_int_conv_VTA_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(SN_path,"pupil_derivative.txt"), deriv_pup_int_conv_SN_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(ACC_path,"pupil_derivative.txt"), deriv_pup_int_conv_ACC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(OCC_path,"pupil_derivative.txt"), deriv_pup_int_conv_OCC_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(BF_path,"pupil_derivative.txt"), deriv_pup_int_conv_BF_roi[0:150], delimiter=',')
            np.savetxt(os.path.join(Pons_path,"pupil_derivative.txt"), deriv_pup_int_conv_Pons_roi[0:150], delimiter=',')
                    
                                        

# 5d. time-to-peak convolution: 
        # convolve with canonical HRF at variable (1-6) TTPs (+pup pushed back 1s)
        pup_int_TR_bins_conv_canonical_P6 = np.convolve(pup_int_TR_bins_center, HRF_canonical, mode='full')
        pup_int_TR_bins_conv_canonical_P5 = np.convolve(pup_int_TR_bins_center, SPM_canonicalHRF_P5, mode='full')
        pup_int_TR_bins_conv_canonical_P4 = np.convolve(pup_int_TR_bins_center, SPM_canonicalHRF_P4, mode='full')
        pup_int_TR_bins_conv_canonical_P3 = np.convolve(pup_int_TR_bins_center, SPM_canonicalHRF_P3, mode='full')
        pup_int_TR_bins_conv_canonical_P2 = np.convolve(pup_int_TR_bins_center, SPM_canonicalHRF_P2, mode='full')
        pup_int_TR_bins_conv_canonical_P1 = np.convolve(pup_int_TR_bins_center, SPM_canonicalHRF_P1, mode='full')
        
# 5d. for the pupil derivative 
        deriv_pup_int_TR_bins_conv_canonical_P6 = np.convolve(deriv_pup_int_TR_center, HRF_canonical, mode='full')
        deriv_pup_int_TR_bins_conv_canonical_P5 = np.convolve(deriv_pup_int_TR_center, SPM_canonicalHRF_P5, mode='full')
        deriv_pup_int_TR_bins_conv_canonical_P4 = np.convolve(deriv_pup_int_TR_center, SPM_canonicalHRF_P4, mode='full')
        deriv_pup_int_TR_bins_conv_canonical_P3 = np.convolve(deriv_pup_int_TR_center, SPM_canonicalHRF_P3, mode='full')
        deriv_pup_int_TR_bins_conv_canonical_P2 = np.convolve(deriv_pup_int_TR_center, SPM_canonicalHRF_P2, mode='full')
        deriv_pup_int_TR_bins_conv_canonical_P1 = np.convolve(deriv_pup_int_TR_center, SPM_canonicalHRF_P1, mode='full')
        
        if save_data:
            save_4_convolved_TTP_dir = os.path.join(raw_datapath, os.listdir(raw_datapath)[subj], 'logfiles-ses-' + day[d], "4_convolved", '4_TTP_analysis')
            if not os.path.exists(save_4_convolved_TTP_dir):
                os.makedirs(save_4_convolved_TTP_dir)
            # save TTP analyses pupil files
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_dilation_canonical_P6.txt"), pup_int_TR_bins_conv_canonical_P6[0:150], delimiter=',')
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_dilation_canonical_P5.txt"), pup_int_TR_bins_conv_canonical_P5[0:150], delimiter=',')
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_dilation_canonical_P4.txt"), pup_int_TR_bins_conv_canonical_P4[0:150], delimiter=',')
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_dilation_canonical_P3.txt"), pup_int_TR_bins_conv_canonical_P3[0:150], delimiter=',')
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_dilation_canonical_P2.txt"), pup_int_TR_bins_conv_canonical_P2[0:150], delimiter=',')
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_dilation_canonical_P1.txt"), pup_int_TR_bins_conv_canonical_P1[0:150], delimiter=',')
    
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_derivative_canonical_P6.txt"), deriv_pup_int_TR_bins_conv_canonical_P6[0:150], delimiter=',')
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_derivative_canonical_P5.txt"), deriv_pup_int_TR_bins_conv_canonical_P5[0:150], delimiter=',')
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_derivative_canonical_P4.txt"), deriv_pup_int_TR_bins_conv_canonical_P4[0:150], delimiter=',')
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_derivative_canonical_P3.txt"), deriv_pup_int_TR_bins_conv_canonical_P3[0:150], delimiter=',') 
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_derivative_canonical_P2.txt"), deriv_pup_int_TR_bins_conv_canonical_P2[0:150], delimiter=',')
            np.savetxt(os.path.join(save_4_convolved_TTP_dir,"pupil_derivative_canonical_P1.txt"), deriv_pup_int_TR_bins_conv_canonical_P1[0:150], delimiter=',')
            
            
   
            
            
            
# save variance of pupil timeseries 
save_all = [save_subj, save_day, save_dat]
save_all=list(map(list, zip(*save_all))) #flip    
save_all_header = ['subject', 'day', 'SD pupil']
save_all.insert(0,save_all_header)
#save_all.tocsv("D:\\NYU_RS_LC\\stats\\R_csv_files\\standard_dev_pupil.csv")




fig, axs = plt.plot()
plt.rc('font', size=22)   
#fig.suptitle('pupil size')
plt.plot.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.plot.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.plot.rc('ytick', labelsize=16)
plt.figure(figsize=(11, 3))
plt.rc('font', size=32) 
#plt.rc('legend',fontsize=32)  
plt.plot(pup_int_TR_bins_center[0:149],'tab:blue', linewidth=1) # downsampled pup, centered, no 1-sec shift \
plt.plot(pup_int_TR_bins_conv_canonical_P1[0:149],'tab:red', linewidth=1)   # downsampled pup, +2s, centered 
plt.plot(pup_int_TR_bins_conv_canonical_murphy[0:149],'tab:grey', linewidth=1)   # convolved with 1s HRF
#plt.plot(sinewave)
plt.ylabel('pupil size')
plt.plot(0, 0, "blue", label="no lag")
plt.plot(0, 0, "red", label = "+2s lag")
plt.plot(0, 0, "green", label = "conv TTP=1")
plt.legend()
#
#start_time = 0
#end_time = 150
#sample_rate = 1
#time = np.arange(start_time, end_time, 1/sample_rate)
#theta = 0
#frequency = 0.1
#amplitude = 500
#sinewave = amplitude * np.sin(2 * np.pi * frequency * time + theta)
#plt.plot(sinewave)

#
#plt.savefig('overlap_plots_new.eps', format='eps')
#
## plots: 
fig, axs = plt.subplots(3)
plt.rc('font', size=22)   
#fig.suptitle('pupil size')
plt.rc('axes', labelsize=22)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=22)    # fontsize of the tick labels
plt.rc('ytick', labelsize=22)
axs[0].plot(l_pup_int,'tab:blue') # raw pup
axs[1].plot(pup_int_TR_bins_center,'tab:red')   # downsampled
axs[2].plot(pup_int_conv_LC_roi[0:150],'tab:green')   # convolved 

#plt.sca(axs[0])
#plt.xticks(np.arange(0, 16000, step=5000))
#plt.yticks([0, 1000, 2000, 3000])
#
#plt.sca(axs[1])
#plt.xticks(np.arange(0, 160, step=50))
#plt.yticks([-1000, 0, 2000])
#
#plt.sca(axs[2])
#plt.xticks(np.arange(0, 160, step=50))
#plt.yticks([-1000, 0, 2000])
#plt.savefig('process_pupsize_plots.eps', format='eps')
#
## plots: 
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
