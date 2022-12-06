function a_fcwml_mri_firstlevel(statsinput)


%% Load paths and files
padi=i_fcwml_infofile(statsinput.SUBJNAME);


%% Get scans [only for 1 pipeline the sessions are seperate]
switch statsinput.which_analysis

    case '0_MakeResiduals'

        if strcmp(statsinput.subanalysis_name,'ses-day1')
            whichscans=fullfile(padi.newBOLDpath{1},{padi.newBOLD{1}.name}');
            which_days=[1];
        elseif strcmp(statsinput.subanalysis_name,'ses-day2')
            whichscans=fullfile(padi.newBOLDpath{2},{padi.newBOLD{2}.name}');
            which_days=[2];
        end

    otherwise

        if statsinput.smoothedata==1
            % For the rest, only some individuals have only one day
            if statsinput.which_session == 9
                whichscans{1}=fullfile(padi.newBOLDpath{1},{padi.snewBOLD{1}.name}');
                whichscans{2}=fullfile(padi.newBOLDpath{2},{padi.snewBOLD{2}.name}');
                which_days=[1 2];
            elseif statsinput.which_session == 1
                whichscans=fullfile(padi.newBOLDpath{1},{padi.snewBOLD{1}.name}');
                which_days=[1];
            elseif statsinput.which_session == 2
                whichscans=fullfile(padi.newBOLDpath{2},{padi.snewBOLD{2}.name}');
                which_days=[2];
            end
        elseif statsinput.smoothedata==0
            % For the rest, only some individuals have only one day
            if statsinput.which_session == 9
                whichscans{1}=fullfile(padi.newBOLDpath{1},{padi.newBOLD{1}.name}');
                whichscans{2}=fullfile(padi.newBOLDpath{2},{padi.newBOLD{2}.name}');
                which_days=[1 2];
            elseif statsinput.which_session == 1
                whichscans=fullfile(padi.newBOLDpath{1},{padi.newBOLD{1}.name}');
                which_days=[1];
            elseif statsinput.which_session == 2
                whichscans=fullfile(padi.newBOLDpath{2},{padi.newBOLD{2}.name}');
                which_days=[2];
            end
        end

end


%% Make nuisance regressors
for c_ses = which_days

    % 1) Get RPs
    R{c_ses}=load(padi.reglog{c_ses}.finalregs);
    RP{c_ses}=R{c_ses}.R(:,27:32);
    RP{c_ses}=RP{c_ses}-RP{c_ses}(1,:); %correction: bc of ME combination, it did not start at 0

    % 2) RETROICOR
    RETROICOR{c_ses}=R{c_ses}.R(:,1:26); %

    % 3) get 4th ventrical signal [from u files]
    inputim.path=padi.BOLDpath{c_ses};
    inputim.ims={padi.BOLD{c_ses}.name}';
    [sigextr{c_ses}, ~] = f_extract_BOLD_data(inputim, padi.r4th_vent ,1);
    sigextr{c_ses}=sigextr{c_ses}'-mean(sigextr{c_ses});

    %     % 4) SCRUBBING [check orginal threshold]
    %     SCRUB{c_ses}=load(padi.reglog{c_ses}.scrub);
    %     spikes = find(~SCRUB{c_ses}.temp_mask);
    %     spikereg_matrix=zeros(length(SCRUB{c_ses}.temp_mask),sum(~SCRUB{c_ses}.temp_mask));
    %     for i=1:numel(spikes)
    %         spikereg_matrix(spikes(i),i)=1;
    %     end
    %     SCRUB{c_ses}=spikereg_matrix;

    % COMBINE
    nuisanceregs{c_ses} = [RP{c_ses} RETROICOR{c_ses} sigextr{c_ses}];

end


%% Make output directory based on (sub)analysis folder
if statsinput.subanalysis == 1
    current_directory=fullfile(padi.stats,'fMRI',...
        statsinput.which_analysis,statsinput.subanalysis_name,statsinput.SUBJNAME);

elseif statsinput.subanalysis == 0
    current_directory=fullfile(padi.stats,'fMRI',...
        statsinput.which_analysis,statsinput.SUBJNAME);
end

%Make directory
if exist(current_directory);rmdir(current_directory,'s');end
mkdir(current_directory)
cd(current_directory)


%% Now do pipeline specific analysis
switch statsinput.which_analysis

    % Now do pipeline specific stuff
    case '0_MakeResiduals' % Does sessions seperate!

        % pupil
        % > not needed

        % save reg file
        if strcmp(statsinput.subanalysis_name,'ses-day1')
            R=nuisanceregs{1};
            save('nuisanceregs_ses-day1.mat','R')
            whichregs=fullfile(current_directory,'nuisanceregs_ses-day1.mat');
        elseif strcmp(statsinput.subanalysis_name,'ses-day2')
            R=nuisanceregs{2};
            save('nuisanceregs_ses-day2.mat','R')
            whichregs=fullfile(current_directory,'nuisanceregs_ses-day2.mat');
        end

    case '1_Murphy_replication'


        for c_ses = which_days

            % pupil
            P1 = load(fullfile(padi.puplog{c_ses},'4_convolved','1_Murphy_replication','pupil_dilation_canonical.txt'));
            P2 = load(fullfile(padi.puplog{c_ses},'4_convolved','1_Murphy_replication','pupil_dilation_dispersion.txt'));
            P3 = load(fullfile(padi.puplog{c_ses},'4_convolved','1_Murphy_replication','pupil_dilation_temporal.txt'));

            % Demean
            P1 = P1 - mean(P1);
            P2 = P2 - mean(P2);
            P3 = P3 - mean(P3);

            % save reg file
            R=[P1 P2 P3 nuisanceregs{c_ses}];
            save(fullfile(current_directory,['allregressors_ses-day',num2str(c_ses),'.mat']),'R')
        end

        % For the rest, only some individuals have only one day
        if statsinput.which_session == 9
            whichregs{1}=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
            whichregs{2}=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        elseif statsinput.which_session == 1
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
        elseif statsinput.which_session == 2
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        end
            

    case '1_Murphy_replication_6RP'

        for c_ses = which_days

            % pupil
            P1 = load(fullfile(padi.puplog{c_ses},'4_convolved','1_Murphy_replication','pupil_dilation_canonical.txt'));
            P2 = load(fullfile(padi.puplog{c_ses},'4_convolved','1_Murphy_replication','pupil_dilation_dispersion.txt'));
            P3 = load(fullfile(padi.puplog{c_ses},'4_convolved','1_Murphy_replication','pupil_dilation_temporal.txt'));

            % Demean
            P1 = P1 - mean(P1);
            P2 = P2 - mean(P2);
            P3 = P3 - mean(P3);

            % save reg file
            R=[P1 P2 P3 RP{c_ses}];
            save(fullfile(current_directory,['allregressors_ses-day',num2str(c_ses),'.mat']),'R')
        end
        
        % For the rest, only some individuals have only one day
        if statsinput.which_session == 9
            whichregs{1}=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
            whichregs{2}=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        elseif statsinput.which_session == 1
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
        elseif statsinput.which_session == 2
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        end

    case '2_Schneider_replication'

        for c_ses = which_days

            % pupil
            P1 = load(fullfile(padi.puplog{c_ses},'4_convolved','2_Schneider_replication','pupil_dilation.txt'));

            % Demean
            P1 = P1 - mean(P1);

            % save reg file
            R=[P1 nuisanceregs{c_ses}];
            save(fullfile(current_directory,['allregressors_ses-day',num2str(c_ses),'.mat']),'R')
        end
        % For the rest, only some individuals have only one day
        if statsinput.which_session == 9
            whichregs{1}=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
            whichregs{2}=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        elseif statsinput.which_session == 1
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
        elseif statsinput.which_session == 2
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        end

    case '3a_rsHRF_estimation'

        % Does not happen in this script!

    case '3b_rsHRF_pupil_dilation'

        for c_ses = which_days

            % pupil
            P1 = load(fullfile(padi.puplog{c_ses},'4_convolved','3_rsHRF',...
                ['pupil_dilation_',statsinput.subanalysis_name,'.txt']));

            % Demean
            P1 = P1 - mean(P1);

            % save reg file
            R=[P1 nuisanceregs{c_ses}];
            save(fullfile(current_directory,['allregressors_ses-day',num2str(c_ses),'.mat']),'R')

        end
        % For the rest, only some individuals have only one day
        if statsinput.which_session == 9
            whichregs{1}=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
            whichregs{2}=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        elseif statsinput.which_session == 1
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
        elseif statsinput.which_session == 2
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        end

    case '3c_rsHRF_pupil_derivative'

        for c_ses = which_days

            % pupil
            P1 = load(fullfile(padi.puplog{c_ses},'4_convolved','3_rsHRF',...
                ['pupil_derivative_',statsinput.subanalysis_name,'.txt']));

            % Demean
            P1 = P1 - mean(P1);

            % save reg file
            R=[P1 nuisanceregs{c_ses}];
            save(fullfile(current_directory,['allregressors_ses-day',num2str(c_ses),'.mat']),'R')

        end
        % For the rest, only some individuals have only one day
        if statsinput.which_session == 9
            whichregs{1}=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
            whichregs{2}=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        elseif statsinput.which_session == 1
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
        elseif statsinput.which_session == 2
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        end

    case '4a_TTP_analysis_estimation'

        % Does not happen in this script!

    case '4b_TTP_analysis_pupil_dilation'

        for c_ses = which_days

            % pupil
            P1 = load(fullfile(padi.puplog{c_ses},'4_convolved','4_TTP_analysis',...
                ['pupil_dilation_canonical_P',statsinput.subanalysis_name,'.txt']));

            % Demean
            P1 = P1 - mean(P1);

            % save reg file
            R=[P1 nuisanceregs{c_ses}];
            save(fullfile(current_directory,['allregressors_ses-day',num2str(c_ses),'.mat']),'R')

        end
        % For the rest, only some individuals have only one day
        if statsinput.which_session == 9
            whichregs{1}=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
            whichregs{2}=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        elseif statsinput.which_session == 1
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
        elseif statsinput.which_session == 2
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        end

    case '4c_TTP_analysis_pupil_derivative'

        for c_ses = which_days

            % pupil
            P1 = load(fullfile(padi.puplog{c_ses},'4_convolved','4_TTP_analysis',...
                ['pupil_derivative_canonical_P',statsinput.subanalysis_name,'.txt']));

            % Demean
            P1 = P1 - mean(P1);

            % save reg file
            R=[P1 nuisanceregs{c_ses}];
            save(fullfile(current_directory,['allregressors_ses-day',num2str(c_ses),'.mat']),'R')
        end
        % For the rest, only some individuals have only one day
        if statsinput.which_session == 9
            whichregs{1}=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
            whichregs{2}=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        elseif statsinput.which_session == 1
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(1),'.mat']);
        elseif statsinput.which_session == 2
            whichregs=fullfile(current_directory,['allregressors_ses-day',num2str(2),'.mat']);
        end

end



%% RUN THE STATS

% Load the correct batch
switch statsinput.which_analysis

    case '0_MakeResiduals'

        % This is done per session
        load batch_file_firstlevel_MakeResiduals

    otherwise

        % For some participants we do not include both sessions
        if statsinput.which_session == 9
            load batch_file_firstlevel
        elseif statsinput.which_session == 1 || statsinput.which_session == 2
            load batch_file_firstlevel_onesession
        end

end

% Change output dir
matlabbatch{1}.spm.stats.fmri_spec.dir = {current_directory};

% Decide which sessions [either both, or just day 1 or day 2]
switch statsinput.which_analysis

    case '0_MakeResiduals'

        %session 1 or 2
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = whichscans;
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {whichregs};

    otherwise

        if statsinput.which_session == 9

            %session 1
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = whichscans{1};
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {whichregs{1}};

            %session 2
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = whichscans{2};
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {whichregs{2}};

        elseif statsinput.which_session == 1 || statsinput.which_session == 2

            %session 1 or 2
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = whichscans;
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {whichregs};

        end

end

% Mask
matlabbatch{1}.spm.stats.fmri_spec.mask = {'/project/3023021.01/TEMP_LINDA/masks/other_masks/SPM_brainmask/brainmask.nii'};

% Change contrasts for the Murphy analysis [rest is already correct in the batch]
switch statsinput.which_analysis

    case '1_Murphy_replication'

        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'canonical';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'dispersion';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 0];
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'temporal';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = ...
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep;

    case '1_Murphy_replication_6RP'

        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'canonical';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'dispersion';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 0];
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'temporal';
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
        matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = ...
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep;

end

% Run batch
spm_jobman('run',matlabbatch);

%% Remove what we dont need

betass=dir(fullfile('beta*.nii'));
for i = 1:numel(betass)
    delete(fullfile(betass(i).name));
end
delete(fullfile('mask.nii'));
delete(fullfile('ResMS.nii'));
delete(fullfile('RPV.nii'));
delete(fullfile('SPM.mat'));


