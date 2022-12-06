clear all; clc;

% Settings
clear all; clc;
homeD='D:\NYU_RS_LC\';
subjpath='D:\NYU_RS_LC\data';
subjlist=dir(fullfile(subjpath,'MRI*'));
addpath('D:\NYU_RS_LC\scripts\1_preproc');
ROI={'DR', 'MR', 'LC', 'VTA', 'SN', 'ACC', 'OCC', 'PONS', 'BF', 'Pons'};

% Run each subject
for c_sess = 1
    for c_subj = 1:70
        for c_roi = 1:numel(ROI)          
            subjlist(c_subj).name
            if c_sess==1 && c_subj == 4
                tSNR='NaN';
            else
                [tSNR] = a_extr_tSNR_perROI(subjlist(c_subj).name, c_sess,ROI{c_roi});
            end
              % SAVE DATA
                    %--------------------------------------------------------------------------

            extract_tsnr(1, 1)=cellstr('subj');
            extract_tsnr(1, c_roi+1)=cellstr(ROI{c_roi});

            extract_tsnr(c_subj+1,1)=cellstr(strcat(num2str(subjlist(c_subj).name)));
            if c_sess==1 && c_subj == 4
                extract_tsnr(c_subj+1,c_roi+1)=cellstr(tSNR);
            else
                extract_tsnr(c_subj+1,c_roi+1)=num2cell(tSNR);
            end
            output_path=fullfile(homeD, 'stats', 'R_csv_files');

            filename=['tSNR_day', num2str(c_sess),'.csv'];
            savefilename=fullfile(output_path,filename);
            cell2csv(savefilename,extract_tsnr);

        end
        
    end
end