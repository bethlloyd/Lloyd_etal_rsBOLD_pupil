clc
clear all

%% reslice 0p89 masks to 2p0 using NN interp

%% Path settings ----------------------------------------------------------
addpath('D:\NYU_RS_LC\scripts\4_ROIs');

LC_mask_folder = 'D:\NYU_RS_LC\ANTs\LC_segmentation_BL\LC_native2template';
LC_mask = dir(fullfile(LC_mask_folder, '2p0_*'));


for mask=1:numel(LC_mask.name)
    
    disp(['checking mask ' LC_mask(mask).name]);

    % reslice 
    load batch_reslice_NN_masks.mat

     %change subject code
    matlabbatch = struct_string_replace(matlabbatch,'2p0_template_space_LC_roi_001.nii',char(LC_mask(mask).name));

    %run batch
    spm_jobman('run',matlabbatch); clear matlabbatch
    
    
end