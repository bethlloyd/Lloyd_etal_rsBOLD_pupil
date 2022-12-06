clear all; clc;
homeD='D:\NYU_RS_LC\';
subjpath='D:\NYU_RS_LC\data';
subjlist=dir(fullfile(subjpath,'MRI*'));
addpath('D:\NYU_RS_LC\scripts');


%% Define looped varaibles ---------------------------------------------------
ROI={'DR', 'MR', 'LC', 'VTA', 'SN', 'ACC', 'OCC', 'PONS', 'BF', 'Pons'};
pup = {'3b_rsHRF_pupil_dilation', '3c_rsHRF_pupil_derivative'};


%for sm = 1
for roi =1:9
    disp(['running..', ROI{roi}]);

    for p = 1:2

        for c_subj=1:numel(subjlist)
            disp(['running..',subjlist(c_subj).name])
            [con_stat] = a_extract_1stlev_stats_rsHRF(subjlist(c_subj).name, ROI{roi}, pup{p});

          % SAVE DATA
            %--------------------------------------------------------------------------
            extract_con_file(1, 1)=cellstr('subj');
            extract_con_file(1, 2)=cellstr('spmT_0001.nii');

            extract_con_file(c_subj+1,1)=cellstr(strcat(num2str(subjlist(c_subj).name)));
            extract_con_file(c_subj+1,2)=num2cell(con_stat);


            stats_dir=fullfile(homeD, 'stats', 'R_csv_files');
            statspath=fullfile(stats_dir, pup{p});       
            
            if ~exist(statspath)
                mkdir(statspath)
            end
            filename=strcat(['con_stat', '_rsHRF_ana_' ROI{roi},'.csv']);
            savefilename=fullfile(statspath,filename);
            cell2csv(savefilename,extract_con_file);

        end
    end 
end
    

