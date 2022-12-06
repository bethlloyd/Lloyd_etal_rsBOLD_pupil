
% EXTRACT LC and PONS SIGNAL INTENSITY FROM INDIV SEGMENTATIONs
%--------------------------------------------------------------------------
% author: BL 2022

% PATH SETTINGS
%--------------------------------------------------------------------------

home='D:\NYU_RS_LC';
addpath('D:\NYU_RS_LC\scripts\0_general');
subjpath=fullfile(home,'data');
subjlist=dir(fullfile(subjpath,'MRI*'));

% Make save folder for intensities
savepath='D:\NYU_RS_LC\stats\R_csv_files\individual_LCs';


% prep stats list
intensity_stats(1,1)=cellstr('subj');
intensity_stats(1,2)=cellstr('LC_intensity');
intensity_stats(1,3)=cellstr('PONS_intensity');
intensity_stats(1,4)=cellstr('LC_cnr_intensity');

% no indiv mask for subj 004, 007, 020 -> exclude these from intensity analysis 
subjexcl=[4, 6, 18]; 
list_total=[1:70];
subjincl=setdiff(list_total,subjexcl);

i = 0
for file=subjincl
    i = i+1
    disp(subjlist(file).name)

    % get subj number
    subjnum=subjlist(file).name(10:end)
    %subjnum = LC_ROI(file).name(10:12);
    
    LC_mask = [subjpath, '\', subjlist(file).name, '\ses-day2\ROI\LC\', 'LC_delin_', subjnum, '.nii']
    if isfile(LC_mask)
        disp('got LC mask')
    else
        disp(['no LC mask for ', subjlist(file).name])
    end
    %get hdr of lc
    lc_hdr=spm_vol(LC_mask);
    
    % load in pons roi 
    pons_roi =[subjpath, '\', subjlist(file).name, '\ses-day2\ROI\Pons\', 'Pons_delin_', subjnum, '.nii']
    
    if isfile(pons_roi)
        disp('got pons mask')
    else
        disp(['no pons mask for ', subjlist(file).name])
    end
    
    pons_hdr=spm_vol(pons_roi);
    
    % load in the fse scan
    FSE_path=dir(['D:\NYU_RS_LC\data\MRI_FCWML', subjnum, '\ses-day2\anat\sub*-fse*']);
    FSEim=fullfile(FSE_path.folder, FSE_path.name,['MRI_FCWML', subjnum, '_0001.nii']);
    fse_hdr=spm_vol(FSEim);

%     %check dimentions
%     if abs(sum(sum(c_hdr.mat-r_hdr.mat)))>0
%         error('ROI and CONTRAST MAP are not in the same space!')
%     end

    %get roi coordinates
    lc_roixyz = f_NYULC_threeDfind(lc_hdr,1);
    pons_roixyz = f_NYULC_threeDfind(pons_hdr,1);
    
    LC_intensity = mean(spm_get_data(fse_hdr,lc_roixyz));
    pons_intensity = mean(spm_get_data(fse_hdr,pons_roixyz));
    LC_cnr = (LC_intensity-pons_intensity)/pons_intensity
    %get number of voxels:
    %num_voxels=numel(spm_get_data(c_hdr,roixyz{c_roi}));

    %get the data from the images based on the ROI
    intensity_stats(i+1,1)=cellstr(strcat(num2str(subjlist(file).name)));
    intensity_stats(i+1,2)=num2cell(LC_intensity);
    intensity_stats(i+1,3)=num2cell(pons_intensity);
    intensity_stats(i+1,4)=num2cell(LC_cnr);
    

end % roi



% SAVE DATA
%--------------------------------------------------------------------------
filename=strcat('LC_intensities_new.csv');
savefilename=fullfile(savepath,filename);
cell2csv(savefilename,intensity_stats);

writecell(intensity_stats, savefilename)