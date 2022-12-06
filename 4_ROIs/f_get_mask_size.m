function vox_no=f_get_mask_size(input_nii)

addpath('D:\NYU_RS_LC\scripts\0_general');

% mask folder 
%mask_folder = ['D:\NYU_RS_LC\masks\masks_template_space_Jun2022', ROI_name];
%mask_file=dir(fullfile(mask_folder, '*.nii'));


%get mask coordinates
roixyz = f_NYULC_threeDfind(input_nii,0.1);

% Calculate mm3
%qMM=size(roixyz,2)*res_EPI
vox_no = size(roixyz,2);
%x = y * z 
disp(vox_no);
%180 = y * 5.832

