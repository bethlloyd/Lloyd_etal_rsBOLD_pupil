function [tSNR] = a_extr_tSNR_perROI(SUBJNAME, sess, roi)
% inputs: 

% SUBJNAME = subject name/code
% session = i.e. day1, day2 (for path definition) 
% len_ims = number of functional images in one session
% roi_mask = path and ROI mask name i.e. ''E:\NYU_RS_LC\masks\Keren mask\rLC_2SD_BINARY_TEMPLATE.nii''
% extrval = set threshold for mask to extract signal 

%% Calculate tSNR 
% this script uses function f_extract_BOLD_data to extract tSNR in
% different brain regions (using ROI masks). 
% Saves a .csv file with a tSNR value (per session, per subject) 


%% path settings 
homeD='D:\NYU_RS_LC\';
addpath('D:\NYU_RS_LC\scripts\0_general');
% roi path
ROI_path = 'D:\NYU_RS_LC\masks\masks_template_space_Jun2022';

% len_ims
len_ims=1;
extrval = 0.1;

%get the tSNR images that are in the subject folder
inputim.path=fullfile(homeD,'stats', 'fMRI', 'tSNR', SUBJNAME, ['ses-day',num2str(sess)]);
BOLDim=dir(fullfile(inputim.path,'tSNR_im.nii'));
inputim.ims={BOLDim.name}';

% check for correct number of func ims
if numel(inputim.ims)~=len_ims
    disp('length of BOLD ims not correct length')
end

%roi path
if strcmp(roi,'LC') == 1
    ROI_folder=dir(fullfile('D:\NYU_RS_LC\data', SUBJNAME, 'ses-day2\ROI\LC', ['r*.nii']));
    ROI_mask=fullfile(ROI_folder.folder, ROI_folder.name);
else
    ROI_folder=dir(fullfile(ROI_path, [roi, '_roi'], ['r*.nii']));
    ROI_mask=fullfile(ROI_folder.folder, ROI_folder.name);
end

% extract signal 
[sigextr, roixyz] = f_extract_BOLD_data(inputim, ROI_mask, extrval);

tSNR=sigextr;




