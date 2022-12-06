function [sigextr, roixyz] = f_extract_BOLD_data(inputim, roitemplate,extrval)

% EXTRACT BOLD SIGNAL WITH ROI MASKS
%--------------------------------------------------------------------------
%
% Input
%   inputim.path = path to the images (beta/raw etc.) you want to extract data from
%   inputim.ims = cell array with the images you want to extract data from;
%   roitemplate = the mask
%   extravl = value of the coordinates in the mask (it will look for 
%   coordinates equal or higher than that value)
%
% Output
%   sigextr = average signal within the mask of the images
%   roixyz = list of coordinates from which extr sig is taken
%
%--------------------------------------------------------------------------


% Get coordinates from the ROI
r_hdr=spm_vol(roitemplate);
roixyz = f_NYULC_threeDfind(r_hdr,extrval);

sigextr=[];
% Loop over BOLD images
for c_bold=1:numel(inputim.ims)
    
        % Get the header
        c_hdr=spm_vol(fullfile(inputim.path,inputim.ims{c_bold}));

%         check dimentions
%         if abs(sum(sum(c_hdr.mat-r_hdr.mat)))>0
%             error('ROI and CONTRAST MAP are not in the same space!')
%         end


        %get the data from the images based on the ROI
        sigextr(c_bold)=nanmean(spm_get_data(c_hdr,roixyz));
       
end

