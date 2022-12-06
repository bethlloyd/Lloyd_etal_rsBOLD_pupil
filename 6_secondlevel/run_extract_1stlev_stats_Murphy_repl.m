clear all; clc;
homeD='D:\NYU_RS_LC\';
subjpath='D:\NYU_RS_LC\data';
subjlist=dir(fullfile(subjpath,'MRI*'));
addpath('D:\NYU_RS_LC\scripts');


%% Define looped varaibles ---------------------------------------------------
ROI={'LC'};
pup = {'1_Murphy_replication'};

HRF1_agg = zeros(1,321)';
HRF2_agg = zeros(1,321)';
HRF3_agg = zeros(1,321)';
%for sm = 1
for roi =1
    disp(['running..', ROI{roi}]);

    for p = 1

        for c_subj=1:numel(subjlist)
            disp(['running..',subjlist(c_subj).name])
            [HRF1, HRF2, HRF3] = a_extract_1stlev_stats_Murphy_repl(subjlist(c_subj).name, ROI{roi}, pup{p});

            HRF1_agg = (HRF1_agg + HRF1);
            HRF2_agg = (HRF2_agg + HRF2);
            HRF3_agg = (HRF3_agg + HRF3);
            
          % SAVE DATA
            %--------------------------------------------------------------------------
%             extract_con_file(1, 1)=cellstr('subj');
%             extract_con_file(1, 2)=cellstr('canonical');
%             extract_con_file(1, 3)=cellstr('temporal');
%             extract_con_file(1, 4)=cellstr('derivative');
%             
%             extract_con_file(c_subj+1,1)=cellstr(strcat(num2str(subjlist(c_subj).name)));
%             extract_con_file(c_subj+1,2)=num2cell(con_stat1);
%             extract_con_file(c_subj+1,3)=num2cell(con_stat2);
%             extract_con_file(c_subj+1,4)=num2cell(con_stat3);
%             
%             stats_dir=fullfile(homeD, 'stats', 'R_csv_files');
%             statspath=fullfile(stats_dir, pup{p});       
%             
%             if ~exist(statspath)
%                 mkdir(statspath)
%             end
%             filename=strcat(['spmT_stat', '_Murphy_repl_' ROI{roi},'.csv']);
%             savefilename=fullfile(statspath,filename);
%             cell2csv(savefilename,extract_con_file);

        end
    end 
end
HRF1_agg=HRF1_agg/numel(subjlist);
HRF2_agg=HRF2_agg/numel(subjlist);
HRF3_agg=HRF3_agg/numel(subjlist);

sum_all = HRF1_agg + HRF2_agg + HRF3_agg;
plot(sum_all)
ylim([-0.004 0.004])