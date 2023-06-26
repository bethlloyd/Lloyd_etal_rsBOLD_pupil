function [dice_c, nvox_overlap ,nvox_rater1, nvox_rater2] = a_calc_dicecoeff(SUBJNAME)

% CALCULATE DICE COEFFICIENT 
%--------------------------------------------------------------------------
% author: BL 2021

% PATH SETTINGS
%--------------------------------------------------------------------------
addpath('D:\NYU_RS_LC\scripts');
addpath('D:\NYU_RS_LC\scripts\0_general');
%get SUBJNAME
if ~exist('SUBJNAME')
    SUBJNAME=char(inputdlg('Which subject?'));
end

num_subj = SUBJNAME(10:12);
%path settings
padi=i_fcwml_infofile(SUBJNAME);


% GET DATA
%--------------------------------------------------------------------------

mask_rater1=padi.rater{1};
mask_rater2=padi.rater{2};

%get the beta map hdr      
c_hdr_rat1=spm_vol(fullfile(mask_rater1.folder,mask_rater1.name));      
c_hdr_rat2=spm_vol(fullfile(mask_rater2.folder,mask_rater2.name));


%check dimentions
if abs(sum(sum(c_hdr_rat1.mat-c_hdr_rat2.mat)))>0
    error('MASKS are not in the same space!')
end

%get roi coordinates
roixyz1 = f_NYULC_threeDfind(fullfile(mask_rater1.folder,mask_rater1.name),.1);
roixyz2 = f_NYULC_threeDfind(fullfile(mask_rater2.folder,mask_rater2.name),.1);

%flip rows and columns 
roixyz1=transpose(roixyz1);
roixyz2=transpose(roixyz2);

%check for equal columns
c_all=[];
for c_r1 = 1:size(roixyz1,1)
    c_all(c_r1)=ismember(roixyz1(c_r1,:),roixyz2,'rows');
end

%log number of vox overlap and number for each rater
n_overlap=sum(c_all);
n_rater1=size(roixyz1, 1);
n_rater2=size(roixyz2, 1);


%caclulcate dice
dice_c = [(2*n_overlap)/(n_rater1+n_rater2)]';
nvox_overlap = [n_overlap]';
nvox_rater1 = [n_rater1]';
nvox_rater2 = [n_rater2]';
%get the data from the images based on the ROI
extract_coeff(1,1)=cellstr(strcat(num2str(SUBJNAME)));
extract_coeff(1,2)=num2cell(dice_c);
     
  
  % SAVE DATA
%--------------------------------------------------------------------------
%remove 2nd row (full of commas) 
extracteddata{c_roi}(2,:) = [];
newroi_name=erase(roinames{c_roi}, '.nii');

% filename=strcat('dice_coeff.csv');
savefilename=fullfile(padi.stats, 'LC_mask',filename);
ell2csv(savefilename,extract_coeff);
% 
% SUBJNAME

end

