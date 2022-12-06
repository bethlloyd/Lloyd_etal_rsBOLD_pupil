function f_check_overlap_masks(SUBJNAME)



% CHECK NO OVERLAP BETWEEN MASKS
%--------------------------------------------------------------------------
% author: BL 2021

% PATH SETTINGS
%--------------------------------------------------------------------------

%get SUBJNAME
if ~exist('SUBJNAME')
    SUBJNAME=char(inputdlg('Which subject?'));
end

%path settings
padi=i_fcwml_infofile(SUBJNAME);

% load in masks
mask_lc = fullfile( padi.LC_roi, 'overlap', 'lcmask_regNat_binary.nii');
mask_4thvent = fullfile(padi.roi, '4th_ventricle', 'r4th_vent.nii');

%get the mask hdr      
hdr_mask_lc=spm_vol(mask_lc);      
hdr_mask_4thvent=spm_vol(mask_4thvent);

%check dimentions
if abs(sum(sum(hdr_mask_lc.mat-hdr_mask_4thvent.mat)))>0
    error('MASKS are not in the same space!')
end

%get roi coordinates
roixyz1 = f_NYULC_threeDfind(mask_lc,.1);
roixyz2 = f_NYULC_threeDfind(mask_4thvent,.1);
disp(['Number of voxels in the mask is: ' num2str(size(roixyz1,2))])
%flip rows and columns 
roixyz1=transpose(roixyz1);
roixyz2=transpose(roixyz2);

%check for equals in  columns
c_all=[];
for c_r1 = 1:size(roixyz1,1)
    c_all(c_r1)=ismember(roixyz1(c_r1,:),roixyz2,'rows');
end

%ALERT for overlapping voxels
if c_all>0
    error('Overlapping voxels!')
    disp(SUBJNAME);
end


