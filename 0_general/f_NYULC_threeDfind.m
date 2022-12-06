function [roixyz] = f_NYULC_threeDfind(scanname,value)
%--------------------------------------------------------------------------
%returns coordinates that above value
%
%ldV2013
%--------------------------------------------------------------------------


%read in ROI
hdr_roitemp=spm_vol(scanname);
data_roitemp=spm_read_vols(hdr_roitemp);

%make ROI search volume (X-Y-X coordinates of 1's in the ROI)
xxsv=zeros(hdr_roitemp.dim);
yysv=zeros(hdr_roitemp.dim);
zzsv=zeros(hdr_roitemp.dim);
for xx=1:hdr_roitemp.dim(1)
    xxsv(xx,:,:)=xx;
end
for yy=1:hdr_roitemp.dim(2)
    yysv(:,yy,:)=yy;
end
for zz=1:hdr_roitemp.dim(3)
    zzsv(:,:,zz)=zz;
end
roixx=xxsv(data_roitemp>=value);
roiyy=yysv(data_roitemp>=value);
roizz=zzsv(data_roitemp>=value);

roixyz=[roixx,roiyy,roizz,ones(size(roixx,1),1)]';