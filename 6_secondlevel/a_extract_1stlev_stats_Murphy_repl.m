function [HRF1, HRF2, HRF3] = a_extract_1stlev_stats_Murphy_repl(SUBJNAME, roi, pup_type)

% EXTRACT BETA SCORES FOR ROI MASKS
%--------------------------------------------------------------------------
% author: BL 2021

% PATH SETTINGS
%--------------------------------------------------------------------------

% SETTINGS
%--------------------------------------------------------------------------
homeD='D:\NYU_RS_LC\';
addpath('D:\NYU_RS_LC\scripts\7_secondlevel');
addpath('D:\NYU_RS_LC\scripts\0_general');
% roi path
ROI_path = 'D:\NYU_RS_LC\masks\masks_template_space_Jun2022';
%get dat dir
statspath=fullfile(homeD, 'stats', 'fMRI', pup_type);

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
con_im1=dir(fullfile(con_dir, 'spmT_0001.nii'));
con_im2=dir(fullfile(con_dir, 'spmT_0002.nii'));
con_im3=dir(fullfile(con_dir, 'spmT_0003.nii'));

%get hdr of roi
r_hdr=spm_vol(ROI_mask);
        
%get the con map hdr
c_hdr1=spm_vol(fullfile(con_im1.folder, con_im1.name));
c_hdr2=spm_vol(fullfile(con_im2.folder, con_im2.name));
c_hdr3=spm_vol(fullfile(con_im3.folder, con_im3.name));

            
%check dimentions
% if abs(sum(sum(c_hdr.mat-r_hdr.mat)))>0
%     error('ROI and CONTRAST MAP are not in the same space!')
% end
            
%get roi coordinates
roixyz = f_NYULC_threeDfind(r_hdr,mask_theshold);
                       
%get the data from the conmap based on the ROI
con_stat1=mean(spm_get_data(c_hdr1,roixyz));
con_stat2=mean(spm_get_data(c_hdr2,roixyz));            
con_stat3=mean(spm_get_data(c_hdr3,roixyz));   


%% make the basis funcs
%---------------------------------------------------------------
dt = 0.1; % in seconds -- this should result in functions with the same temporal
% resolution as the main dataset
% then several excerts from spm_get_bf.m
% canonical hemodynamic response function
%---------------------------------------------------------------
p = [6 16 1 1 6 0 32];
bf = spm_hrf(dt,p);

dp = 1;
p(6) = p(6) + dp;
D = (bf(:,1) - spm_hrf(dt,p))/dp;
bf = [bf D(:)];
p(6) = p(6) - dp;
% add dispersion derivative
%--------------------------------------------------------
dp = 0.01;
p(3) = p(3) + dp;
D = (bf(:,1) - spm_hrf(dt,p))/dp;
bf = [bf D(:)];


% multply constat by the bf
HRF1 = con_stat1*bf(:,1);
HRF2 = con_stat2*bf(:,2);
HRF3 = con_stat3*bf(:,3);



