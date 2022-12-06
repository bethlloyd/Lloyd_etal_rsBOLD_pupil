# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 09:36:45 2018

@author: linda
"""
import numpy as np
import re
import matplotlib.pyplot as plt

#libraries
class EyeLinkRead:
    
    """Do analysis on eyelink data"""
    
    def __init__(self):
        
        self.eyelink_sF=1000
        self.sF=250
        
        self.rawdat = []
        self.pupdat = []
        
        # Interpolation settings
        self.win_ave=100 # how many samples used for averaging
        self.win_lim=200 # minimum difference between eye blink events
        
        #turn vis on off
        self.vis=0
        

    def do_readfile(self, filename):
        
        """Load file"""
        
        if self.pupdat==[]:
        
            # read file
            f = open(filename, 'r')
            for c,line in enumerate(f.readlines()):
                self.rawdat.append(line)
            
            # remove lines [script crashes when these lines are in there]
            for rmstr in ["EFIX","ESACC","SFIX","SSACC","SBLINK","EBLINK","END"]:
                self.rawdat = [line for line in self.rawdat if not rmstr in line]
            
            # take the pupil diameter
            dat=[ [], [], [], [] ]
            self.rawevt=[]
            
            #get recording line so everything  before will be removed
            rcd_line=[cnt for cnt, line in enumerate(self.rawdat) if "SAMPLES" in line]

            #get starttime [used to subtract from event times]
            #self.starttime=int(self.rawdat[rcd_line[0]+1].split()[0])
            
            # remove those lines from rawdata
            self.rawdat = self.rawdat[rcd_line[-1]+1:]
            
            # get events and make vector
            for c,line in enumerate(self.rawdat):
                
                #if c>rcd_line[-1]: # just throws away lines before calibration [but might be file specific]

                if "MSG" in line:
                    self.rawevt.append(line)
                else:
                    spt_line=line.split('\t')
                    for cc,dd in enumerate(spt_line):
                        if cc<4:
                            dat[cc].append(dd.strip('  '))
            
            #get pupil dilation
            self.pupdat = [int(float(x)) for x in dat[3]]
            
        else:
            print("Cannot read file twice!")

    def do_cutoutdata(self,starttime,endtime):
        
        """ Cut out a piece of the data to have the desired timecourse length
        for example to match the length of your fMRI data
        """
        
        #get starttime [used to subtract from event times]
        self.starttime=starttime
        
        self.pupdat = self.pupdat[starttime:endtime]

    def get_eyeblinks(self):
        
        """Get eye blinks for interpolation"""
                       
        # make empty
        self.interpol_evt_str=[]
        self.interpol_evt_end=[]

        # Get null events
        interpol_vec=[int(x==0) for x in self.pupdat] 
                      
        # Do not do this for the beginning and end
        interpol_vec[0:self.win_lim]=[0]*self.win_lim
        interpol_vec[-self.win_lim+1:-1]=[0]*self.win_lim
        
        # Get start and end of each interpolation window
        for c, num in enumerate(np.diff(interpol_vec)):
            if num == 1:
                self.interpol_evt_str.append(c+1)
            elif num == -1:
                self.interpol_evt_end.append(c)
        
        # Remove last start value if there is more than end
        if len(self.interpol_evt_str)>len(self.interpol_evt_end):
            self.interpol_evt_str.pop()
        
        #remove "eyeblinks" occure to close in time
        for c, num in enumerate(self.interpol_evt_str): 
            if c<len(self.interpol_evt_str)-1:
                if (self.interpol_evt_str[c+1]-self.interpol_evt_end[c])<self.win_lim:
                    self.interpol_evt_str[c+1]=0
                    self.interpol_evt_end[c]=0
               
        #remove zero's        
        self.interpol_evt_str=[num for num in self.interpol_evt_str if num!=0]
        self.interpol_evt_end=[num for num in self.interpol_evt_end if num!=0]
        
        # OLD methods, did not seem to work well...
#       # Check "eyeblinks" occure to close in time
#       dif_str=[0]
#       dif_str.extend([int(x) for x in np.diff(self.interpol_evt_str)<self.win_lim])
#       self.interpol_evt_str=[num for c,num in enumerate(self.interpol_evt_str) if dif_str[c]==0]
#       dif_end=dif_str[1:]+[0]
#       self.interpol_evt_end=[num for c,num in enumerate(self.interpol_evt_end) if dif_end[c]==0]
        
        
            
        
    def get_largedeviations(self):
        
        """Maybe there are other deviations you may want to interpolate"""
        
        pass
        
    
    def do_interpol(self):
        
        """Interpolate the eye blinks"""
        
        #copy
        self.int_pupdat=self.pupdat[:]
        
        #treat first differently when smaller that the window
        if self.interpol_evt_str[0]<self.win_ave:
            
            #temp window length till end of data
            win_ave_end=len(self.pupdat)-self.interpol_evt_str[-1]-1
            
             # Define interpolation value [average of -X window]
            str_val=self.pupdat[self.interpol_evt_str[0]-1]
            end_val=self.pupdat[self.interpol_evt_end[0]+self.win_ave]
            
            # Define the gap that needs to be filled [half the start window to half the end window]
            gap_val=(self.interpol_evt_end[0]+self.win_ave)-(self.interpol_evt_str[0]-1)
            int_val=(end_val-str_val)/gap_val
            
            # interpolate
            for c_sam in range((self.interpol_evt_str[0]-1),(self.interpol_evt_end[0]+self.win_ave)+1):
                self.int_pupdat[c_sam]=self.int_pupdat[c_sam-1]+int_val
            
            #remove from the list
            self.interpol_evt_str.pop(0)
            self.interpol_evt_end.pop(0)
        
        #treat last differently when smaller that the window
        if len(self.pupdat)-self.interpol_evt_end[-1]<self.win_ave:
            
            #temp window length till end of data
            win_ave_end=len(self.pupdat)-self.interpol_evt_end[-1]-1
            
            # Define interpolation value [average of -X window]
            str_val=self.pupdat[self.interpol_evt_str[-1]-self.win_ave]
            end_val=self.pupdat[self.interpol_evt_end[-1]+win_ave_end]
            
            # Define the gap that needs to be filled [half the start window to half the end window]
            gap_val=(self.interpol_evt_end[-1]+win_ave_end)-(self.interpol_evt_str[-1]-self.win_ave)
            int_val=(end_val-str_val)/gap_val
        
            # interpolate
            for c_sam in range((self.interpol_evt_str[-1]-self.win_ave),(self.interpol_evt_end[-1]+win_ave_end)+1):
                self.int_pupdat[c_sam]=self.int_pupdat[c_sam-1]+int_val
                
            self.interpol_evt_str.pop()
            self.interpol_evt_end.pop()
        
        # Loop over interpolation start values
        for c_evt, n_evt in enumerate(self.interpol_evt_str):
            
            # Define interpolation value [average of -X window]
            str_val=self.pupdat[self.interpol_evt_str[c_evt]-self.win_ave]
            end_val=self.pupdat[self.interpol_evt_end[c_evt]+self.win_ave]
            
            # Define the gap that needs to be filled [half the start window to half the end window]
            gap_val=(self.interpol_evt_end[c_evt]+self.win_ave)-(self.interpol_evt_str[c_evt]-self.win_ave)
            int_val=(end_val-str_val)/gap_val
        
            # interpolate
            for c_sam in range((self.interpol_evt_str[c_evt]-self.win_ave),(self.interpol_evt_end[c_evt]+self.win_ave)+1):
                self.int_pupdat[c_sam]=self.int_pupdat[c_sam-1]+int_val
    
    def get_events(self,eventnames):
        
        """Get event markers"""
        
        self.eventtimes=[]
        
        if type(eventnames)!=list:
            print("input should be a list!")
            
        else:
        
            #loop over events
            for gettrig in eventnames:
                
                #get raw values from rawdat
                rawvals=[int(re.sub("[^0-9]", "",line)) for line in self.rawdat if gettrig in line]
                
                #subtract starttime so it is in sync with pupdat
                downsF=int(self.eyelink_sF/self.sF)
                self.eventtimes.append([int((num-self.starttime)/downsF) for num in rawvals])
                
    def plot_events(self,n_evt,plot_win):
        
        """Plot the events"""
        
        # Plot each trial seperate
        fig, axs = plt.subplots(int(np.ceil(len(self.eventtimes[n_evt])/3)),3, figsize=(15, 12), facecolor='w', edgecolor='k')
        fig.subplots_adjust(hspace = .5, wspace=.5)
        
        # Loop over trials
        for ax, d in zip(axs.ravel(), self.eventtimes[n_evt][0:len(self.eventtimes[n_evt])]):
            
            #plot the trial
            ax.plot(self.int_pupdat[d-(self.sF*np.abs(plot_win[0])):d+(self.sF*plot_win[1])])  

        plt.show()
        
    def do_man_interpol(self,startvals,endvals):
        
        """You may want to manually interpolate some pieces of data"""
        
        # Loop over interpolation start values
        for c_evt, n_evt in enumerate(startvals):
            
            # Define interpolation value [average of -X window]
            str_val=self.int_pupdat[startvals[c_evt]-self.win_ave]
            end_val=self.int_pupdat[endvals[c_evt]+self.win_ave]
            
            # Define the gap that needs to be filled [half the start window to half the end window]
            gap_val=(endvals[c_evt]+self.win_ave)-(startvals[c_evt]-self.win_ave)
            int_val=(end_val-str_val)/gap_val
        
            # interpolate
            for c_sam in range((startvals[c_evt]-self.win_ave),(endvals[c_evt]+self.win_ave)+1):
                self.int_pupdat[c_sam]=self.int_pupdat[c_sam-1]+int_val
    
    
    def get_pupresp(self,n_evt,trial_win,resp_win,bsl_win):
        
        """Get the pupil responses and trial data for later group analysis"""
        
        self.trialdata=np.array(list(range(1,((trial_win[1]-trial_win[0])*self.sF)+1)))
        self.pupresponses=[]
  
        for c_tr, d in enumerate(self.eventtimes[n_evt]):
            
            if d > self.sF*np.abs(bsl_win[0]):
            
                #baseline to peak collection
                bl=np.mean(self.int_pupdat[d-(self.sF*np.abs(bsl_win[0])):d-(self.sF*np.abs(bsl_win[1]))])
                pk=np.mean(self.int_pupdat[d+(self.sF*resp_win[0]):d+int(self.sF*resp_win[1])]/bl)
                self.pupresponses.append(pk)
            
                #collect inter trial for later averaging
#                EE.trialdata=np.append(EE.trialdata, np.array([EE.int_pupdat[d-(EE.sF*np.abs(trial_win[0])):
#                    d+(EE.sF*trial_win[1])]])/bl, axis=0)
                self.trialdata=np.vstack([self.trialdata,np.array([self.int_pupdat[d-(self.sF*np.abs(trial_win[0])):d+(self.sF*trial_win[1])]])/bl])
    
            else:
                self.pupresponses.append(np.nan)
                self.trialdata=np.vstack([self.trialdata,[np.nan]*(trial_win[1]-trial_win[0])*self.sF])

            
            

        
        
#
        
        
            
                
            
            
        
        
#        
#    #get triggers
#stimtrig=[]
#for gettrig in ["CS_ONSET"]:
#    starttime=[int(re.sub("[^0-9]", "",line)) for line in rawdat if "START_TIME" in line]
#    endtime=[int(re.sub("[^0-9]", "",line)) for line in rawdat if "END_TIME" in line]
#    CS_onsets = [int(re.sub("[^0-9]", "",line)) for line in rawdat if gettrig in line]
#CS_onsets=[int((line - starttime[0])/4) for line in CS_onsets]
#print("This file has " + str(len(CS_onsets)) + " CS events")
#
#
#
#
#
#
#
#def interpol2(pupdat):
#    
#    # Interpolation settings
#    win_ave=50 # how many samples used for averaging
#    win_lim=200 # minimum difference between eye blink events
#    int_std=6
#    
#    # Create new variable which will be the interpolated data
#    int_pupdat=pupdat[:]
#    
#    # Define interpolation [if val_to_int == 0 then only blinks]
#    val_to_int=np.max(int_pupdat)-(np.std(pupdat)*int_std)
#    interpol_vec=[int(x<val_to_int) for x in pupdat] 
#    
#    # Get start and end of each interpolation window
#    interpol_evt_str=[]
#    interpol_evt_end=[]
#    for c, num in enumerate(np.diff(interpol_vec)):
#        if num == 1:
#            interpol_evt_str.append(c+1)
#        elif num == -1:
#            interpol_evt_end.append(c)
#            
#            
#    # Check "eyeblinks" occure to close in time
#    dif_str=[0]
#    dif_str.extend([int(x) for x in np.diff(interpol_evt_str)<win_lim])
#    interpol_evt_str=[num for c,num in enumerate(interpol_evt_str) if dif_str[c]==0]
#    dif_end=dif_str[1:]+[0]
#    interpol_evt_end=[num for c,num in enumerate(interpol_evt_end) if dif_end[c]==0]
#    
#    # Loop over interpolation start values
#    for c_evt, n_evt in enumerate(interpol_evt_str):
#        
#        # Define interpolation value [average of -X window]
#        str_val=pupdat[interpol_evt_str[c_evt]-win_ave]
#        end_val=pupdat[interpol_evt_end[c_evt]+win_ave]
#        
#        # Define the gap that needs to be filled [half the start window to half the end window]
#        gap_val=(interpol_evt_end[c_evt]+win_ave)-(interpol_evt_str[c_evt]-win_ave)
#        int_val=(end_val-str_val)/gap_val
#    
#        for c_sam in range((interpol_evt_str[c_evt]-win_ave),(interpol_evt_end[c_evt]+win_ave)+1):
#            int_pupdat[c_sam]=int_pupdat[c_sam-1]+int_val
#    
#    plt.plot(pupdat)
#    plt.plot(int_pupdat)
#
#
#    return int_pupdat