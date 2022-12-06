function a_rsHRF_est(SUBJNAME)

% Paths and files
padi = i_fcwml_infofile(SUBJNAME);

% current directory
current_directory = fullfile('/project/3023021.01/TEMP_LINDA/stats/fMRI/3a_rsHRF_estimation/', SUBJNAME);
mkdir(current_directory)
cd(current_directory)

% adjust the regressor file to add together day 1 and day 2
regs_D1 = load(padi.reglog{1}.finalregs).R;
regs_D2 = load(padi.reglog{2}.finalregs).R;
R = [regs_D2];% regs_D2];
save(fullfile(current_directory,['allregressors_concat.mat']),'R')

% Load
load 1_canonical_batchfile.mat

% Change input
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.images = [fullfile(padi.newBOLDpath{2},{padi.newBOLD{2}.name}')]%;...
    %fullfile(padi.newBOLDpath{2},{padi.newBOLD{2}.name}')];

matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.Denoising.generic{1}.multi_reg = {fullfile(current_directory,['allregressors_concat.mat'])};

groupims={matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.genericROI{1}.images{2:end}}';
matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.genericROI{1}.images = [{padi.indvLC} matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.genericROI{1}.images{2:end}]';

matlabbatch{1}.spm.tools.rsHRF.ROI_rsHRF.outdir = {current_directory};

% Run batch
spm_jobman('run',matlabbatch);