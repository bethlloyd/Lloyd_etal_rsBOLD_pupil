function [Cxy,F_mscohere, Pxy,F_cpsd]=f_fcwml_crossSpectrum(SUBJNAME, pup_type, roi, session)


% get the denoised BOLD (concatonated then seperate)
% load in the pupil vectors
% perform cross-spectrum for each session

%% SETTINGS

%path settings
addpath('D:\NYU_RS_LC\scripts');
addpath('D:\NYU_RS_LC\scripts\0_general');
padi=i_fcwml_infofile(SUBJNAME);

% settings
formatSpec = '%f';
extrval=0.01;
%rng default
Window_length = 10;
no_overlap = 3;
Fs = 0.5;
%freqrange = 'twosided';


%% remove excluded subjects
subjexcl_d1={'MRI_FCWML004', 'MRI_FCWML016', 'MRI_FCWML025', 'MRI_FCWML043', 'MRI_FCWML055',... 
                'MRI_FCWML119', 'MRI_FCWML160', 'MRI_FCWML219'};
subjexcl_d2={'MRI_FCWML020', 'MRI_FCWML022', 'MRI_FCWML028', 'MRI_FCWML033', 'MRI_FCWML036', 'MRI_FCWML064'};
which_day = {subjexcl_d1, subjexcl_d2};         

% residual data folder 
RES_dat = 'D:\NYU_RS_LC\stats\fMRI\0_MakeResiduals';
% roi path
ROI_path = 'D:\NYU_RS_LC\masks\masks_template_space_Jun2022';

size_len = 31;

if ismember(SUBJNAME, which_day{session})  
    
%     Cxy= NaN(1,len_F,'single');
%     F_mscohere= NaN(1,len_F,'single');
%     Pxy= NaN(1,len_F,'single');
%     F_cpsd= NaN(1,len_F,'single');
    Cxy = NaN(1,size_len,'single');
    F_mscohere= NaN(1,size_len,'single');
    
    Pxy = NaN(1,size_len,'single');
    F_cpsd= NaN(1,size_len,'single');
    
    
else
   % Load in RESIDUAL data
    inputim.path=fullfile(RES_dat,padi.sessions{session},SUBJNAME);

    BOLDim=dir(fullfile(inputim.path,['Res_*']));
    inputim.ims={BOLDim.name}';
    
    %roi path
    if strcmp(roi,'LC') == 1
        ROI_folder=dir(fullfile('D:\NYU_RS_LC\data', SUBJNAME, 'ses-day2\ROI\LC', ['r_2p0*.nii']));
        ROI_mask=fullfile(ROI_folder.folder, ROI_folder.name);
    else
        ROI_folder=dir(fullfile(ROI_path, [roi, '_roi'], ['r_2p0*.nii']));
        ROI_mask=fullfile(ROI_folder.folder, ROI_folder.name);
    end
    % extract signal
    [sigextr_roi, roixyz] = f_extract_BOLD_data(inputim, ROI_mask,extrval);
    %sigextr_roi=sigextr_roi'-mean(sigextr_roi); % demean the signal

    
      % assign the BOLD signal 
    if session == 1
        sess = 'ses-day1';
        %BOLD_sig = HRF_dat.data(1:150,roi);
    elseif session == 2
        sess = 'ses-day2';
        %BOLD_sig = HRF_dat.data(151:300,roi);
    end
    
      %load in the pupil vectors
    %load in the pupil vectors
    pup_dir=fullfile('D:\NYU_RS_LC\data', SUBJNAME, ['logfiles-',sess], '3_downsampled');
    
    if pup_type==1% Load in pupil regressor # 2 - Diameter
    % get correct folder
        pup_filename=['pupil_dilation_1s_shifted.txt'];
        pupdat=fullfile(pup_dir, pup_filename);
        
        % try to open
        [fid, errormsg] = fopen(pupdat, 'r+');
        if ~isempty(errormsg)
            warning('failed to open %s due to %s', pupdat, errormsg);
        else
             % open datafile
            pup_int = fscanf(fid,formatSpec);
            fclose(fid);
        end
        
    elseif pup_type==2 % Load in pupil regressor # 2 - Derivative
        % get correct folder
        deriv_pup_filename=['pupil_derivative_1s_shifted.txt'];
        deriv_pupdat=fullfile(pup_dir, deriv_pup_filename);
        
        % try to open
        [fid, errormsg] = fopen(deriv_pupdat, 'r+');
        if ~isempty(errormsg)
            warning('failed to open %s due to %s', deriv_pupdat, errormsg);
        else
             % open datafile
            pup_int = fscanf(fid,formatSpec);
            fclose(fid);
        end

    end
    
    pup_int=pup_int'-mean(pup_int); % demean the signal & zscore both timeseries
    pup_int= zscore(pup_int');
    sigextr_roi = zscore(sigextr_roi');
    
    [Cxy,F_mscohere] = mscohere(pup_int, sigextr_roi, Window_length, no_overlap, 60, Fs);

%     plot(F_mscohere,Cxy)
%     title('Magnitude-Squared Coherence')
%     xlabel('Frequency (Hz)')
%     grid
    %plot(Cxy)
%     % Cross Spectrum 
    [Pxy,F_cpsd] = cpsd(pup_int,sigextr_roi, Window_length, no_overlap, 60, Fs);
    %Pxy(Cxy < 0.01) = 0;
    
    % normalise (demean) the timcourse 
     %Cxy_demeaned=Cxy - mean(Cxy);
    
%     plot(F_cpsd,Pxy)
%     title('Cross Spectrum Phase')
%     xlabel('Frequency (Hz)')
%     ylabel('Lag (\times\pi rad)')
%     grid

%     plot(F_cpsd,angle(Pxy_orig)/pi)
%     title('Cross Spectrum Phase')
%     xlabel('Frequency (Hz)')
%     ylabel('Lag (\times\pi rad)')
%     grid
end



