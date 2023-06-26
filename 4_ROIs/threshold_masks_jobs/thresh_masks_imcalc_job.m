%-----------------------------------------------------------------------
% Job saved on 23-May-2022 17:03:41 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.util.imcalc.input = {'F:\NYU_RS_LC\MNI_space_masks_FINAL\OCC_roi\AAL\1p8_MNI_calcarine_L+R.nii,1'};
matlabbatch{1}.spm.util.imcalc.output = '30_thresh_1p8_MNI_calcarine_L+R.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {'F:\NYU_RS_LC\MNI_space_masks_FINAL\OCC_roi\AAL'};
matlabbatch{1}.spm.util.imcalc.expression = 'i1>0.3';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
