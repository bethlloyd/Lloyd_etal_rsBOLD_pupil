clear all; clc;
homeD='D:\NYU_RS_LC\';
subjpath='D:\NYU_RS_LC\data';
subjlist=dir(fullfile(subjpath,'MRI*'));
addpath('D:\NYU_RS_LC\scripts');


%% Define looped varaibles ---------------------------------------------------
ROI={'DR', 'MR', 'LC', 'VTA', 'SN', 'ACC', 'OCC', 'PONS', 'BF', 'Pons'};
models={'1', '2', '3', '4', '5', '6'};
%smooth = {'smoothed', 'unsmoothed'};
pup = {'4b_TTP_analysis_pupil_dilation', '4c_TTP_analysis_pupil_derivative'};

for mod = 1:numel(models)

    %for sm = 1
    for roi =1:9
        disp(['running..', ROI{roi}]);

        for p = 1:2

            for c_subj=1:numel(subjlist)
                disp(['running..',subjlist(c_subj).name])
                [con_stat] = a_extract_1stlev_stats_TTP(subjlist(c_subj).name, ROI{roi}, models{mod}, pup{p});
         
              % SAVE DATA
                %--------------------------------------------------------------------------
                extract_con_file(1, 1)=cellstr('subj');
                extract_con_file(1, 2)=cellstr('spmT_0001.nii');

                extract_con_file(c_subj+1,1)=cellstr(strcat(num2str(subjlist(c_subj).name)));
                extract_con_file(c_subj+1,2)=num2cell(con_stat);


                stats_dir=fullfile(homeD, 'stats', 'R_csv_files');
                statspath=fullfile(stats_dir, pup{p});       

                filename=strcat(['con_stat', '_rsHRF_P', models{mod}, '_' ROI{roi},'.csv']);
                savefilename=fullfile(statspath,filename);
                cell2csv(savefilename,extract_con_file);

            end
        end 
    end
    
end

