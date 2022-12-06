clc
clear all

%% Path settings ----------------------------------------------------------
addpath('D:\NYU_RS_LC\scripts\4_ROIs');

LC_mask_folder = 'D:\NYU_RS_LC\ANTs\LC_segmentation_BL\LC_native2template';
LC_mask = dir(fullfile(LC_mask_folder, '0p89_*'));

% Collect
hdr = {'Subjname', 'nVox','Vox_mm3'};

% Loop
for mask=1:numel(LC_mask)
    
    disp(['checking mask ' LC_mask(mask).name])
    
    load batch_reslice_NN_masks
    
    %change subject code
    matlabbatch = struct_string_replace(matlabbatch,'0p89_template_space_LC_roi_001.nii',char(LC_mask(mask).name));
    matlabbatch = struct_string_replace(matlabbatch,'MRI_FCWML001',['MRI_FCWML' LC_mask(mask).name(end-6:end-4)]);
    
    %run batch
    spm_jobman('run',matlabbatch);
    
    %print
    newfile = fullfile(matlabbatch{2}.spm.util.imcalc.outdir{1},matlabbatch{2}.spm.util.imcalc.output);
    [roixyz] = f_NYULC_threeDfind2(newfile,0,2);
    
    alldat(mask,1) = num2cell(size(roixyz,2));
    alldat(mask,2) = num2cell(size(roixyz,2)*2*2*2);
    
    clear matlabbatch
    
    
end

% Save csv file
newdat = [{LC_mask.name}' alldat];
writecell([hdr; newdat],'LC_nVoxelSize.csv')


%% GROUP LC

for mask=1:numel(LC_mask)
    
    disp(['checking mask ' LC_mask(mask).name])
    
    load batch_groupLC_NN
    
    %change subject code
    matlabbatch = struct_string_replace(matlabbatch,'0p89_template_space_LC_roi_001.nii',char(LC_mask(mask).name));
    
    matlabbatch{1}.spm.util.imcalc.output = char(fullfile(matlabbatch{1}.spm.util.imcalc.outdir,['bin_' LC_mask(mask).name]));
    
    %run batch
    spm_jobman('run',matlabbatch);
    
    clear matlabbatch
    
end

%% Create group mask for num 4,7,20

scanname = 'agg_0p89_LC_ROIs_69p.nii';

[roixyz] = f_NYULC_threeDfind2(scanname,1,2);

hdr = spm_vol(scanname);
        
dat = spm_get_data(hdr,roixyz);

dat(dat > (mean(dat) + (std(dat)*2)));

% Use outcome variable in imcalc to create 2SD mask
