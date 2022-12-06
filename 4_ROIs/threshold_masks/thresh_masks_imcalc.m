% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'E:\NYU_RS_LC\scripts\4_ROIs\threshold_masks\thresh_masks_imcalc_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
