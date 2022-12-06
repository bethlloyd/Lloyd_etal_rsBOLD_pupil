
%--------------------------------------------------------------------------
%
% perform 1st level for FCWML
%
%BL 2021
%--------------------------------------------------------------------------

statsdir='D:\NYU_RS_LC\stats\fMRI';
whichana={'1_Murphy_replication', '1_Murphy_replication_6RP'};

for c_ana = 1:numel(whichana)
    
    current_directory=fullfile(statsdir,whichana{c_ana});

    cd(current_directory)
    mkdir(fullfile(current_directory,'groupstats','F-test'))
    
    %load batch
    load batch_second_level_Murphy_F_test
    
    matlabbatch = struct_string_replace(matlabbatch,'1_Murphy_replication',whichana{c_ana});

%     %get confiles
%     conpath=fullfile(current_directory,'groupstats','T_Pos');
%     confiles=cellstr(spm_select('List',conpath,['^.*\.nii']));
%     
%     %change input
%     matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile(current_directory,'groupstats','F-test')};
%     matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = ...
%         strcat(conpath,filesep,confiles);
%     %change subject code

    %run job
    spm_jobman('run',matlabbatch); clear matlabbatch

end


