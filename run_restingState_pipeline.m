clear all;
clc;
%--------------------------------------------------------------------------

%GENERAL SETTINGS
%--------------------------------------------------------------------------

% Path defenition
project.mainpath='/project/3023021.01/TEMP_LINDA/';
project.datapath=fullfile(project.mainpath,'data');

% Add scripts path
addpath(genpath(fullfile(project.mainpath,'scripts')));

% Add fieldtrip to use parallel computing
addpath('/home/common/matlab/fieldtrip/qsub')

% Which steps
do_ants=0;
do_smooth=0;
do_frstlevel=0;
do_rsHRF=1;
do_cluster=0;

%% Subject settings

% Subject info
subjlist = {'MRI_FCWML001', 'MRI_FCWML002', 'MRI_FCWML003', 'MRI_FCWML004', 'MRI_FCWML006',... %1-5
    'MRI_FCWML007', 'MRI_FCWML008', 'MRI_FCWML009', 'MRI_FCWML010', 'MRI_FCWML011',... %6-10
    'MRI_FCWML012', 'MRI_FCWML013', 'MRI_FCWML014', 'MRI_FCWML015', 'MRI_FCWML016',... %11-15
    'MRI_FCWML017', 'MRI_FCWML018', 'MRI_FCWML020', 'MRI_FCWML022', 'MRI_FCWML023',... %16-20
    'MRI_FCWML025', 'MRI_FCWML026', 'MRI_FCWML027', 'MRI_FCWML028', 'MRI_FCWML029',... %21-25
    'MRI_FCWML030', 'MRI_FCWML031', 'MRI_FCWML032', 'MRI_FCWML033', 'MRI_FCWML035',... %26-30
    'MRI_FCWML036', 'MRI_FCWML037', 'MRI_FCWML038', 'MRI_FCWML040', 'MRI_FCWML041',... %31-35
    'MRI_FCWML042', 'MRI_FCWML043', 'MRI_FCWML046', 'MRI_FCWML047',...                 %36-39
    'MRI_FCWML048', 'MRI_FCWML049', 'MRI_FCWML051', 'MRI_FCWML052', 'MRI_FCWML053',... %40-44
    'MRI_FCWML054', 'MRI_FCWML055', 'MRI_FCWML056', 'MRI_FCWML057', 'MRI_FCWML061',... %45-49
    'MRI_FCWML063', 'MRI_FCWML064', 'MRI_FCWML104', 'MRI_FCWML107', 'MRI_FCWML108',... %50-54
    'MRI_FCWML109', 'MRI_FCWML119', 'MRI_FCWML121', 'MRI_FCWML128', 'MRI_FCWML129',... %55-59
    'MRI_FCWML134', 'MRI_FCWML150', 'MRI_FCWML158', 'MRI_FCWML159',...                 %60-63
    'MRI_FCWML160', 'MRI_FCWML162', 'MRI_FCWML205', 'MRI_FCWML219', 'MRI_FCWML255',... %64-68
    'MRI_FCWML320', 'MRI_FCWML328'};
subjincl = [4];

subjlist_which_session = [9 9 9 2 9, ...
    9 9 9 9 9, ...
    9 9 9 9 2, ...
    9 9 1 1 9, ...
    2 9 9 1 9, ...
    9 9 9 1 9, ...
    1 9 9 9 9, ...
    9 2 9 9, ...
    9 9 9 9 9, ...
    9 2 9 9 9, ...
    9 1 9 9 9, ...
    9 2 9 9 9, ...
    9 9 9 9, ...
    2 9 9 2 9, ...
    9 9]; %9=both, 1=day1,2=day2

% Stats info
analysis_list = {...
    '0_MakeResiduals',...
    '1_Murphy_replication',...
    '1_Murphy_replication_6RP',...
    '2_Schneider_replication',...
    '3b_rsHRF_pupil_dilation',...
    '3c_rsHRF_pupil_derivative',...
    '4b_TTP_analysis_pupil_dilation',...
    '4c_TTP_analysis_pupil_derivative'};
analincl=[2];
subanalysis_bool = [1 0 0 0 1 1 1 1];
smoothedata_bool = [0 1 1 1 0 0 0 0];
subanalysis_list = {...
    {'ses-day1', 'ses-day2'},...
    [],[],[],...
    {'LC', 'VTA', 'SN', 'DR', 'MR', 'BF', 'ACC', 'OCC', 'Pons'},...
    {'LC', 'VTA', 'SN', 'DR', 'MR', 'BF', 'ACC', 'OCC', 'Pons'},...
    {'1', '2', '3', '4', '5', '6'},...
    {'1', '2', '3', '4', '5', '6'}};

%% 1) ANTS
%--------------------------------------------------------------------------
if do_ants == 1

    for c_ses = 1:2

        % First level
        for c_subj = subjincl

            %go to seperate folder to dump job files
            cd /project/3023021.01/TEMP_LINDA/scripts/TORQUEJOBS

            %submit job to cluster
            strname=[subjlist{c_subj} '_day' num2str(c_ses) '_ants'];
            qsubfeval(@a_dotransforms,...
                subjlist{c_subj},...       %input number 1
                c_ses,...                   %input number 2
                'memreq',1000*(1024^2),...
                'timreq',500*60,...
                'batchid',strname,...
                'diary','always')

        end

    end

end

%% 2) SMOOTH 6mm
%--------------------------------------------------------------------------
if do_smooth == 1

    for c_ses = 1:2

        % First level
        for c_subj = subjincl

            %go to seperate folder to dump job files
            cd /project/3023021.01/TEMP_LINDA/scripts/TORQUEJOBS

            %submit job to cluster
            strname=[subjlist{c_subj} '_day' num2str(c_ses) '_smooth'];
            if do_cluster == 1
                qsubfeval(@a_smooth,...
                    subjlist{c_subj},...       %input number 1
                    c_ses,...                   %input number 2
                    'memreq',1000*(1024^2),...
                    'timreq',500*60,...
                    'batchid',strname,...
                    'diary','always')
            else
                disp(strname)
                a_smooth(subjlist{c_subj},c_ses);
            end

        end

    end

end

%% 3) First level subanalincl
%--------------------------------------------------------------------------

if do_frstlevel == 1

    % Loop over analysis [there are 10]
    for c_ana = analincl

        % First level
        for c_subj = subjincl

            % Does it have a subanalysis? [ses-day1, or LC, or 1TTP?]
            if subanalysis_bool(c_ana) == 1

                % Yes, then loop over subanalysis
                for c_subana = 1:numel(subanalysis_list{c_ana})

                    strname=[subjlist{c_subj} '_' analysis_list{c_ana} '_' subanalysis_list{c_ana}{c_subana}];

                    statsinput.SUBJNAME = subjlist{c_subj};
                    statsinput.which_session = subjlist_which_session(c_subj);
                    statsinput.which_analysis = analysis_list{c_ana};
                    statsinput.subanalysis =subanalysis_bool(c_ana);
                    statsinput.subanalysis_name = subanalysis_list{c_ana}{c_subana};
                    statsinput.smoothedata=smoothedata_bool(c_ana);

                    %go to seperate folder to dump job files
                    cd /project/3023021.01/TEMP_LINDA/scripts/TORQUEJOBS

                    %submit job to cluster
                    if do_cluster == 1
                        qsubfeval(@a_fcwml_mri_firstlevel,...
                            statsinput,...       %input number 1
                            'memreq',5000*(1024^2),...
                            'timreq',500*60,...
                            'batchid',strname,...
                            'diary','always')
                    else
                        disp(strname)
                        a_fcwml_mri_firstlevel(statsinput)
                    end

                end


            elseif subanalysis_bool(c_ana) == 0

                strname=[subjlist{c_subj} '_' analysis_list{c_ana}];

                % No then do regular analysis
                statsinput.SUBJNAME = subjlist{c_subj};
                statsinput.which_session = subjlist_which_session(c_subj);
                statsinput.which_analysis = analysis_list{c_ana};
                statsinput.subanalysis = subanalysis_bool(c_ana);
                statsinput.subanalysis_name = [];
                statsinput.smoothedata=smoothedata_bool(c_ana);

                %go to seperate folder to dump job files
                cd /project/3023021.01/TEMP_LINDA/scripts/TORQUEJOBS

                %submit job to cluster
                if do_cluster == 1
                    qsubfeval(@a_fcwml_mri_firstlevel,...
                        statsinput,...       %input number 1
                        'memreq',5000*(1024^2),...
                        'timreq',500*60,...
                        'batchid',strname,...
                        'diary','always')
                else
                    disp(strname)
                    a_fcwml_mri_firstlevel(statsinput)
                end

            end
        end
    end
end

%% rsHRF
%--------------------------------------------------------------------------
if do_rsHRF == 1

    % First level
    for c_subj = subjincl

        %go to seperate folder to dump job files
        cd /project/3023021.01/TEMP_LINDA/scripts/TORQUEJOBS

        %submit job to cluste
        strname=[subjlist{c_subj} '_rshrf'];
        if do_cluster == 1
        qsubfeval(@a_rsHRF_est,...
            subjlist{c_subj},...       %input number 1
            'memreq',1000*(1024^2),...
            'timreq',500*60,...
            'batchid',strname,...
            'diary','always')
        else
            disp(strname)
            a_rsHRF_est(subjlist{c_subj})
        end

    end

    

end


