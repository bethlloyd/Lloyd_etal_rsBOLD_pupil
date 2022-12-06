clear all; clc;

% Settings
%subjpath='C:\Users\lloydb\surfdrive\ExperimentData\NYU_RS_LC\data';
subjpath='D:\NYU_RS_LC\data';
subjlist=dir(fullfile(subjpath,'MRI*'));


subjexcl={'020', '007', '004'};


% Run each subject
for c_subj = 1:numel(subjlist)
    subjlist(c_subj).name;
    num_subj = subjlist(c_subj).name(10:12);
    
    if ismember(num_subj, subjexcl);
        disp(['skipping', num_subj]);
        continue
    end
    
    [dice_c, nvox_overlap ,nvox_rater1, nvox_rater2] = a_calc_dicecoeff(subjlist(c_subj).name);
    %disp(nvox_overlap)

  % SAVE DATA
    %--------------------------------------------------------------------------
    extract_coeff_file(c_subj,1)=cellstr(strcat(num2str(subjlist(c_subj).name)));
    extract_coeff_file(c_subj,2)=num2cell(dice_c);
    extract_coeff_file(c_subj,3)=num2cell(nvox_overlap);
    extract_coeff_file(c_subj,4)=num2cell(nvox_rater1);
    extract_coeff_file(c_subj,5)=num2cell(nvox_rater2);
    
    
    %extract_coeff_file(2,:) = [];
    filename=strcat('LC_dice_coeff.csv');
    savefilename=fullfile('D:\NYU_RS_LC\stats\R_csv_files\individual_LCs',filename);
    cell2csv(savefilename,extract_coeff_file);

    
    
end
    
