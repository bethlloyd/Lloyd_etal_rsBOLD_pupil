# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 11:00:12 2022

@author: lloydb
"""

#Library
from nilearn import plotting
import matplotlib.pyplot as plt
import os

projectdir = 'D:\\NYU_RS_LC'
statsdir = os.path.join(projectdir,'stats','fMRI','unsmoothed_TTP', '4c_TTP_analysis_pupil_derivative')
outfilepath = os.path.join(projectdir,'stats','figures', 'TTP_analysis')
os.chdir(outfilepath)

#Settings
stat_thrs=2.6
anat_img = os.path.join(projectdir,'grouptemplates','step1_TemplateMultivariateSyN', '0p89_straight_T_template0.nii')

modnames = ['1', '2', '3', '4', '5', '6']
tmapnames = ['respmT_0001.nii','respmT_0001.nii','respmT_0001.nii', 'respmT_0001.nii', 'respmT_0001.nii', 'respmT_0001.nii']

fig, (axes1,axes2,axes3,axes4,axes5, axes6) = plt.subplots(6,1, figsize=(12,16),facecolor=(0, 0, 0))
for (mod, tmap, axes) in zip(modnames, tmapnames,[axes1,axes2,axes3,axes4,axes5, axes6]):
    
    #Get tmap
    tmap_filename = os.path.join(statsdir,mod,'groupstats',
                                 'one-samplettest',tmap)
    
     #Plot
    display = plotting.plot_stat_map(tmap_filename, 
                    #title=f'time-to-peak: {mod} s',
                    threshold=stat_thrs,
                    bg_img=anat_img,
                    draw_cross=False,
                    vmax=stat_thrs+3,
                    symmetric_cbar='True',
                    black_bg=True,
                    dim=0,
                    display_mode="z",
                    annotate=False,
                    axes=axes,
                    cut_coords=[0.10-12,-4,4,12])
    #display.annotate(size=6)
    #display.title(f'time-to-peak: {mod} s', x=0.01, y=1.1, size=8, color='white', bgcolor='black')

plt.show()
display.savefig(f'TTPana_stats_pup_deriv.jpeg')


# plot Murphy and schneider 
Murphy_dir = os.path.join(projectdir, 'stats', 'fMRI', '1_Murphy_replication')
stat_thrs_Ftest=5.626
mni_anat = 'D:\\NYU_RS_LC\\grouptemplates\\Temp2MNI2009b\\mni_icbm152_t1_tal_nlin_asym_09b_hires.nii'
#Get Fmap
Fmap_filename = os.path.join(Murphy_dir,'groupstats',
                             'F-test','mni_spmF_0006.nii')
display = plotting.plot_stat_map(Fmap_filename, 
                    #title=f'time-to-peak: {mod} s',
                    threshold=stat_thrs_Ftest,
                    bg_img=mni_anat,
                    draw_cross=False,
                    vmax=stat_thrs_Ftest,
                    symmetric_cbar='True',
                    black_bg=True,
                    dim=0,
                    display_mode="z",
                    annotate=False,
                    cut_coords=[-15,-10,-5,0,5,10,15])
display.annotate(size=6)
display.title(f'Murphy replication', x=0.01, y=1.1, size=8, color='white', bgcolor='black')
