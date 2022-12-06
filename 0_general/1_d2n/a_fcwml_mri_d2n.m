function a_fcwml_mri_d2n(SUBJNAME)
%--------------------------------------------------------------------------
%
% dicom to niftii convertion using SPM [MRI_FCWML]
%
%LdV2019
%--------------------------------------------------------------------------


%get SUBJNAME
if ~exist('SUBJNAME')
    SUBJNAME=char(inputdlg('Which subject?'));
end


%path settings
padi=i_fcwml_infofile(SUBJNAME);
dicompath=padi.rawdata;
datapath=padi.data;

%when downloading from the CBI calander the "_" gets removed
SUBJNAMEcbi=['sub-' erase(SUBJNAME,'_')];


% Do not change anything below [only bugs]
%--------------------------------------------------------------------------

for c_ses = 1:numel(padi.sessions)

    % Loop over dicom folders [in 'raw']
    for c_type = 1:numel(padi.datatypes)
        
        %get dir paths
        sourcedir=fullfile(dicompath,SUBJNAME,padi.sessions{c_ses},'sourcedata',...
            SUBJNAMEcbi,padi.sessions{c_ses},padi.datatypes{c_type});
        newdir=fullfile(datapath,SUBJNAME,padi.sessions{c_ses},padi.datatypes{c_type});
        
        %make dir
        mkdir(newdir)

        %get files
        datafiles=dir(fullfile(sourcedir,'*.tgz'));

        for c_files = 1:numel(datafiles)
            
            if ~contains(datafiles(c_files).name,'sbref') && ~contains(datafiles(c_files).name,'rest')

                % copy to avoid permission errors
                copyfile(fullfile(sourcedir,datafiles(c_files).name),...
                    fullfile(newdir,datafiles(c_files).name));

                % cd [if not if gives an error with untar function...]
                cd(newdir)

                % untar files so they are dcm files 
                % NOTES: 'untar' on matlab2017 and MAC gives errors
                % therefore if it does not work with try/catch use unix 
                % 'tar' function
                try untar(fullfile(newdir,datafiles(c_files).name))
                catch
                    eval(['!tar -zxf ' fullfile(newdir,datafiles(c_files).name)])
                end

                %path
                outputdir=fullfile(newdir,datafiles(c_files).name(1:end-10));

                %get dicom images
                dicomfiles=dir([outputdir,filesep,'*.dcm']);

                %load SPM job [does the dicom to niftii convertion]
                load f_d2n

                %change input of the SPM BATCH
                matlabbatch{1}.spm.util.import.dicom.data = ...
                    cellstr(strcat(fullfile(outputdir,filesep,cellstr(char(dicomfiles.name)))));
                matlabbatch{1}.spm.util.import.dicom.outdir = ...
                    cellstr(outputdir);

                %run the SPM BATCH
                spm_jobman('run',matlabbatch);

                %get niftiis
                newfiles=dir(fullfile(outputdir,'*.nii'));

                %loop over niftiis and change name
                for c_f=1:numel(newfiles)

                    oldfile=fullfile(outputdir,newfiles(c_f).name);

                    if c_f<10
                        newfile=fullfile(outputdir,...
                            [char(SUBJNAME) '_000' num2str(c_f) '.nii']);
                    elseif c_f<100
                        newfile=fullfile(outputdir,...
                            [char(SUBJNAME) '_00' num2str(c_f) '.nii']);
                    elseif c_f<1000
                        newfile=fullfile(outputdir,...
                            [char(SUBJNAME) '_0' num2str(c_f) '.nii']);
                    else
                        newfile=fullfile(outputdir,...
                            [char(SUBJNAME) '_' num2str(c_f) '.nii']);
                    end

                    %by moving the file the name is changed
                    movefile(oldfile,newfile)

                end %c_files

                % delete dcm
                for c_dcm=1:numel(dicomfiles)
                    delete(fullfile(outputdir,dicomfiles(c_dcm).name))
                end

                % delete tar file
                f=dir('*.tgz');
                for c_f=1:numel(f)
                    delete(fullfile(newdir,f(c_f).name))
                end

                clear matlabbatch
                
            end %do not for ref scans

        end %c_files

    end %c_type
    
end %sessions

