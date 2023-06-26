function f_fcwml_calc_tsnr(SUBJNAME,sess)

stats='D:\NYU_RS_LC\stats\fMRI\tSNR_denoised';
padi=i_fcwml_infofile(SUBJNAME);
%func path
funcs = padi.BOLD_aff{sess};
funcs_path = padi.BOLDpath{sess};

%resid path
RES_dat = 'D:\NYU_RS_LC\stats\fMRI\0_MakeResiduals';
inputim.path=fullfile(RES_dat,padi.sessions{sess},SUBJNAME);
BOLDim=dir(fullfile(inputim.path,['Res_*']));
%create output dir
warning off
mkdir(fullfile(stats,SUBJNAME,['ses-day', num2str(sess)]));
warning on
       
%make cell array
P=fullfile(inputim.path,'/',{BOLDim.name})';

%run tsnr scripts
Q = f_calc_tsnr(P);

%move image
copyfile(Q,...
    fullfile(stats,SUBJNAME,['ses-day', num2str(sess)],'tSNR_im.nii'))        



