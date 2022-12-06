
%% This script collects all the HRF vectors 

%% Path settings ----------------------------------------------------------
%Clus_data='/data/lloydb/data';
homeD='D:\NYU_RS_LC\';
stats=fullfile(homeD,'stats');
datpath=fullfile(stats, 'fMRI', '3a_rsHRF_estimation');
subjlist=dir(fullfile(datpath,'MRI*'));

addpath('D:\NYU_RS_LC\scripts');
addpath('D:\NYU_RS_LC\scripts/0_general');

%% Define subject numbers -------------------------------------------------
%subj_list = 70;


%% Define ROIs
roi = {'LC', 'VTA', 'SN', 'DR', 'MR', 'BF', 'ACC', 'OCC', 'Pons'};

% save headers
%--------------------------------------------------------------------------
extract_HRF(1,1)=cellstr('subj');
extract_HRF(1,2)=cellstr('time');
for c_roi = 1:numel(roi)
    extract_HRF(1,2+(c_roi))=cellstr(roi{c_roi});
end

% save headers = event number
%--------------------------------------------------------------------------
event_num(1,1)=cellstr('subj');
for c_roi = 1:numel(roi)
    event_num(1,1+(c_roi))=cellstr(roi{c_roi});
end

% save headers = TTP of the BOLD response
TTP(1,1)=cellstr('subj');
for c_roi = 1:numel(roi)
    TTP(1,1+(c_roi))=cellstr(roi{c_roi});
end


%% GATHER HRF VECTORS -----------------------------------------------------
% Loop over subjects
for c_subj = 1:numel(subjlist)
    subjlist(c_subj).name
    
    for c_roi = 1:numel(roi)
        
        %% load in HRF scruct file
        hrf_filename = ['rsHRF_est_aff_u', subjlist(c_subj).name,'_0006_hrf.mat'];
        hrf_vec = load(fullfile(datpath, subjlist(c_subj).name,hrf_filename));
        % get the event numbers for each ROI 
        event_num_all=hrf_vec.event_number(c_roi);
        
        % get the x coording (timing)
        len_sec = 32;
        div_factor = numel(hrf_vec.hrfa(:,1))/len_sec;
        time_step = 1/div_factor;
        x = 0:time_step:len_sec;
        x(1)=[];
        
        % get HRF response
        y_all=hrf_vec.hrfa(:,c_roi);

        % unit-normalize the response
        y_hrf_norm = y_all - y_all(1);  % subtract value at t1
        y_hrf_norm=y_hrf_norm/max(y_hrf_norm);   % divide timepoints by maximum value
        y_hrf = y_hrf_norm';


        % get the HRF for each roi 
        if c_subj==1
            extract_HRF(c_subj+1:numel(x)+1,1)=cellstr(strcat(num2str(subjlist(c_subj).name)));
            extract_HRF(c_subj+1:numel(x)+1,2)=num2cell(x);
            extract_HRF(c_subj+1:numel(x)+1,2+(c_roi))= num2cell(y_hrf);
            
            %for the TTP collection
            subj_HRF = extract_HRF(c_subj+1:numel(x)+1,:);
        else
            %get the correct start and end rows
            start_row = numel(x)*(c_subj-1)+2;
            end_row = numel(x)*(c_subj-1)+2+numel(x)-1;

            %save the HRF vectors to struct
            extract_HRF(start_row:end_row,1)=cellstr(strcat(num2str(subjlist(c_subj).name)));
            extract_HRF(start_row:end_row,2)=num2cell(x);
            extract_HRF(start_row:end_row,2+(c_roi))=num2cell(y_hrf);
            
            %for the TTP collection
            subj_HRF = extract_HRF(start_row:end_row,:);
        end
        
         %save to struct
        event_num(c_subj+1, 1)=cellstr(strcat(num2str(subjlist(c_subj).name)));
        event_num(c_subj+1, 1+(c_roi))=num2cell(event_num_all);

        %% save event onsets
        disp(['saving .. ', subjlist(c_subj).name]);

    end %roi
    
   %% get time to peak (based on max BOLD response)
    [max_values,idx]=max(cell2mat(subj_HRF(2:end,3:11)));
    out=[cell2mat(subj_HRF(idx',2)) max_values'];
%     
    for c_roi = 1:numel(roi)
        TTP(c_subj+1, 1)=cellstr(strcat(num2str(subjlist(c_subj).name)));
        TTP(c_subj+1,1+(c_roi))= num2cell(out(c_roi,1));
    end %roi
  
    
    

end %subject


%% SAVE FILES -------------------------------------------------------------

%save extracted HRF
filename=strcat('rsHRF_extracted_HRFs_both_days.csv');
savefilename=fullfile(stats, 'R_csv_files',filename);
cell2csv(savefilename,extract_HRF);
% % 

% save event numbers
filename=strcat('rsHRF_roi_event_num_both_days.csv');
savefilename=fullfile(stats, 'R_csv_files', filename);
cell2csv(savefilename,event_num);

%save time to peak values
filename=strcat('rsHRF_roi_TTP_both_days.csv');
savefilename=fullfile(stats, 'R_csv_files',filename);
cell2csv(savefilename,TTP);

