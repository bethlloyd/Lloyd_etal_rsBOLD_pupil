clear all; clc;

% Path settings
home='F:\NYU_RS_LC\';
homeE='E:\NYU_RS_LC\';
addpath('D:\NYU_RS_LC\scripts');

subjpath='D:\NYU_RS_LC\data';
subjlist=dir(fullfile(subjpath,'MRI*'));
list_total=[1:70];

ROIs = {'DR', 'MR', 'VTA', 'ACC', 'LC', 'SN', 'OCC', 'BF', 'Pons'};
day = {'ses-day1', 'ses-day2'};
pup = {'pup_size', 'pup_deriv'};


for d = 1:2

    for roi=1:9 %1:numel(ROIs)

        for p = 1:2

            for c_subj = list_total

                disp(['running...', subjlist(c_subj).name, ROIs{roi}, day{d}]);
                % SAVE DATA
                %--------------------------------------------------------------------------
                extract_crosscorr(1, 1)=cellstr('subj');
                extract_crosscorr(1, 2)=cellstr('lag');
                extract_crosscorr(1, 3)=cellstr('CC');

                lag = 9;
                start_row = lag*(c_subj-1)+2;
                end_row =lag*(c_subj-1)+2+lag-1;
                lags = [-4;-3; -2; -1; 0; 1; 2; 3; 4];
                [corr]=f_fcwml_crossCorrelations(subjlist(c_subj).name, ROIs{roi}, p, d);

                if c_subj == 1                
                    extract_crosscorr(c_subj+1:numel(lags)+1,1)=cellstr(strcat(num2str(subjlist(c_subj).name)));
                    extract_crosscorr(c_subj+1:numel(lags)+1,2)=num2cell(lags);
                    extract_crosscorr(c_subj+1:numel(lags)+1,3)=num2cell(corr);

                else

                    extract_crosscorr(start_row:end_row,1)=cellstr(strcat(num2str(subjlist(c_subj).name)));
                    extract_crosscorr(start_row:end_row,2)=num2cell(lags);
                    extract_crosscorr(start_row:end_row,3)=num2cell(corr);
                end

            end 
            
            output_path=fullfile('D:\NYU_RS_LC\stats\R_csv_files\crossCorr',day{d});
            % make output dir if none
            if ~exist(output_path, 'dir')
                mkdir(output_path);
            end

            filename=strcat(['1s_pupshift_check_XCorr_stat_', ROIs{roi}, '_', pup{p} '.csv']);
            savefilename=fullfile(output_path,filename);
            cell2csv(savefilename,extract_crosscorr);


        end
       


        %end 
    end
end
                
                
                
                
                
                