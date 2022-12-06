function [con_stat] = a_extract_1stlev_stats_TTP(SUBJNAME, roi, TTP_mod, pup_type)

% EXTRACT BETA SCORES FOR ROI MASKS
%--------------------------------------------------------------------------
% author: BL 2021

% PATH SETTINGS
%--------------------------------------------------------------------------
homeD='D:\NYU_RS_LC\';
addpath('D:\NYU_RS_LC\scripts\7_secondlevel');
addpath('D:\NYU_RS_LC\scripts\0_general');
% roi path
ROI_path = 'D:\NYU_RS_LC\masks\masks_template_space_Jun2022';

%get dat dir
statspath=fullfile(homeD, 'stats', 'fMRI','unsmoothed_TTP' ,pup_type, TTP_mod);

%roi path
if strcmp(roi,'LC') == 1
    ROI_folder=dir(fullfile('D:\NYU_RS_LC\data', SUBJNAME, 'ses-day2\ROI\LC', ['r*.nii']));
    ROI_mask=fullfile(ROI_folder.folder, ROI_folder.name);
else
    ROI_folder=dir(fullfile(ROI_path, [roi, '_roi'], ['r*.nii']));
    ROI_mask=fullfile(ROI_folder.folder, ROI_folder.name);
end

mask_theshold=0.1;

%% Get con stat data 
con_dir=fullfile(statspath, SUBJNAME);
con_im=dir(fullfile(con_dir, 'spmT_0001.nii'));

%get hdr of roi
r_hdr=spm_vol(ROI_mask);
        
%get the con map hdr
c_hdr=spm_vol(fullfile(con_im.folder, con_im.name));
            
%check dimentions
% if abs(sum(sum(c_hdr.mat-r_hdr.mat)))>0
%     error('ROI and CONTRAST MAP are not in the same space!')
% end
            
%get roi coordinates
roixyz = f_NYULC_threeDfind(r_hdr,mask_theshold);
                       
%get the data from the conmap based on the ROI
con_stat=mean(spm_get_data(c_hdr,roixyz));
            

