function [r p] = a_extract_BS_signal_and_correlate(SUBJNAME,sess)

% ADD REGRESSORS TO MAKE FINAL .MAT FILE
%--------------------------------------------------------------------------
% author: BL 2021

% PATH SETTINGS
%--------------------------------------------------------------------------
addpath('D:\NYU_RS_LC\scripts');
addpath('D:\NYU_RS_LC\scripts\0_general');

%path settings
padi=i_fcwml_infofile(SUBJNAME);

%define the BS ROIs 
ROI={'LC', 'VTA', 'SN', 'DR', 'MR', 'BF'};

% residual data folder 
RES_dat = 'D:\NYU_RS_LC\stats\fMRI\0_MakeResiduals';

% GET DATA
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------       
% Load in RESIDUAL data
inputim.path=fullfile(RES_dat,padi.sessions{sess},SUBJNAME);

BOLDim=dir(fullfile(inputim.path,['Res_*']));

ROI_path = 'D:\NYU_RS_LC\masks\masks_template_space_Jun2022';
inputim.ims={BOLDim.name}';

extrval=0.01;
for roi = 1:numel(ROI)
    
    %roi path
    if strcmp(ROI{roi},'LC') == 1
        ROI_folder=dir(fullfile('D:\NYU_RS_LC\data', SUBJNAME, 'ses-day2\ROI\LC', ['r_2p0*.nii']));
        ROI_mask=fullfile(ROI_folder.folder, ROI_folder.name);
    else
        ROI_folder=dir(fullfile(ROI_path, [ROI{roi}, '_roi'], ['r_2p0*.nii']));
        ROI_mask=fullfile(ROI_folder.folder, ROI_folder.name);
    end

    % extract signal
    [sigextr_roi, roixyz] = f_extract_BOLD_data(inputim, ROI_mask,extrval);
    sigextr_roi=sigextr_roi'-mean(sigextr_roi); % demean the signal
    ROI_sig_struct(1:150,roi)=sigextr_roi';
  
end
% pons ROI
ROI_folder=dir(fullfile(ROI_path,  'PONS_roi', ['r_2p0_*.nii']));
ROI_mask=fullfile(ROI_folder.folder, ROI_folder.name);
% extract signal
[sigextr_roi, roixyz] = f_extract_BOLD_data(inputim, ROI_mask,extrval);
sigextr_roi_pons=sigextr_roi'-mean(sigextr_roi); % demean the signal


%correlated vectors 
[r p] = partialcorr(ROI_sig_struct, sigextr_roi_pons);





