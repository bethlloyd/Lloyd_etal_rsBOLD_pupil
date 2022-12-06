
%--------------------------------------------------------------------------
%
% perform 1st level for FCWML
%
%BL 2021
%--------------------------------------------------------------------------

statsdir='/project/3023021.01/TEMP_LINDA/stats/fMRI_smoothed_3mm';
whichana='4c_TTP_analysis_pupil_derivative';
subdirs={'1', '2', '3', '4', '5', '6'};

for c_dir = 1:numel(subdirs)
    
    current_directory=fullfile(statsdir,whichana,subdirs{c_dir});

    cd(current_directory)
    mkdir(fullfile(current_directory,'groupstats','one-samplettest'))
    
    %load batch
    load batch_second_level
    
    %get confiles
    conpath=fullfile(current_directory,'groupstats','T_Pos');
    confiles=cellstr(spm_select('List',conpath,['^.*\.nii']));
    
    %change input
    matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile(current_directory,'groupstats','one-samplettest')};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = ...
        strcat(conpath,filesep,confiles);
    
    %run job
    spm_jobman('run',matlabbatch); clear matlabbatch

end


