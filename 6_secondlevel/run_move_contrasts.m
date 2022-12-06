% Add the pupil regressor to the "final_regressors.mat" file which included
% the realigment parameters and RETROICORplus regressors

clear all; clc;

% Settings
subjlist ={'MRI_FCWML001', 'MRI_FCWML002', 'MRI_FCWML003', 'MRI_FCWML004', 'MRI_FCWML006',... %1-5
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

statsdir='D:\NYU_RS_LC\stats\fMRI';

% need to just run this now for rs-HRF ana
% Stats info
analysis_list = {...
    '1_Murphy_replication',...
    '1_Murphy_replication_6RP',...
    '2_Schneider_replication',...
    '3b_rsHRF_pupil_dilation',...
    '3c_rsHRF_pupil_derivative',...
    '4b_TTP_analysis_pupil_dilation',...
    '4c_TTP_analysis_pupil_derivative'};
subanalysis_bool = [0 0 0 1 1 1 1];
subanalysis_list = {...
    [],[],[],...
    {'LC', 'VTA', 'SN', 'DR', 'MR', 'BF', 'ACC', 'OCC', 'Pons'},...
    {'LC', 'VTA', 'SN', 'DR', 'MR', 'BF', 'ACC', 'OCC', 'Pons'},...
    {'1', '2', '3', '4', '5', '6'},...
    {'1', '2', '3', '4', '5', '6'}};


% subject settings
run_subs = [1:70];

for c_ana = 1:numel(analysis_list)
    
    if subanalysis_bool(c_ana) == 1
        
        for sub_ana = 1:numel(subanalysis_list{c_ana})
            current_directory=fullfile(statsdir, 'unsmoothed_TTP', analysis_list{c_ana}, subanalysis_list{c_ana}{sub_ana});
            mkdir(fullfile(current_directory,'groupstats','T_Pos'))
            
            
        end
    else
        
        current_directory=fullfile(statsdir,analysis_list{c_ana});
        mkdir(fullfile(current_directory,'groupstats','T_Pos'))
        
        % Loop over subjects
        for c_subj = run_subs

            %get number of contrasts:
            conpath=fullfile(current_directory, subjlist{c_subj});
            confile=fullfile(conpath,'con_0001.nii');

            if isfile(confile)
                % move file
                movefile(...
                    confile,...
                    fullfile(current_directory,'groupstats','T_Pos',...
                    [subjlist{c_subj} '_con_0001.nii']))
            end
        end
    end
end
