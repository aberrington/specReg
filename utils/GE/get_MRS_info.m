function [SV_info]=get_MRS_info(SV_info)



%% get data
if ~nargin
    [SV_file,SV_path]=uigetfile('*','get MRS P-file','');
    SV_info=read_MR_headers([SV_path,SV_file],'all');
    disp('######################')
    disp([SV_file,':'])
end

% some info

disp('######################')
disp(['TE = ',num2str(SV_info.image.te/1e3),'ms'])
disp(['TR = ',num2str(SV_info.image.tr/1e3),'ms'])
disp(['Transmitter gain = ',num2str(SV_info.rdb_hdr.ps_mps_tg)])
disp(['receiver gain 1 = ',num2str(SV_info.rdb_hdr.ps_mps_r1)])
disp(['receiver gain 2 = ',num2str(SV_info.rdb_hdr.ps_mps_r2)])
disp(['Volume = ',num2str((SV_info.rdb_hdr.roilenx*SV_info.rdb_hdr.roileny*SV_info.rdb_hdr.roilenz)/1e3),'ml'])
disp(['Number of spectra = ',num2str(SV_info.rdb_hdr.user4)])




%% get MRS orient info
%those are RAS coordinates and give centre of voxel?!
tlhc_LPS(1)=-SV_info.image.tlhc_R;
tlhc_LPS(2)=-SV_info.image.tlhc_A;
tlhc_LPS(3)=SV_info.image.tlhc_S;

trhc_LPS(1)=-SV_info.image.trhc_R;
trhc_LPS(2)=-SV_info.image.trhc_A;
trhc_LPS(3)=SV_info.image.trhc_S;

brhc_LPS(1)=-SV_info.image.brhc_R;
brhc_LPS(2)=-SV_info.image.brhc_A;
brhc_LPS(3)=SV_info.image.brhc_S;

e1_SVS_n=trhc_LPS-tlhc_LPS;
e1_SVS_n=e1_SVS_n./norm(e1_SVS_n);
e2_SVS_n=brhc_LPS-trhc_LPS;
% e2_SVS_n=trhc_LPS-brhc_LPS;
e2_SVS_n=e2_SVS_n./norm(e2_SVS_n);
e3_SVS_n=-1*cross(e1_SVS_n,e2_SVS_n);




disp('orientation of the MRS volume:')
[orientation_SVS]=check_orient(e3_SVS_n);

% if orientation_SVS~=3
%     warndlg('MRS volume is not axial. Please double check correct orientation');
% end

LPS_SVS_centre(1)=-SV_info.image.user11;%L
LPS_SVS_centre(2)=-SV_info.image.user12;%P
LPS_SVS_centre(3)=SV_info.image.user13;%S

% NEED TO CHECK IF THAT IS TRUE FOR NON AXIAL IMAGING VOLUMES
if orientation_SVS==3 %axial
    
    e1_FOV_mm_SVS=SV_info.image.user8;%LR
    e2_FOV_mm_SVS=SV_info.image.user9;%AP
    e3_FOV_mm_SVS=SV_info.image.user10;%CC
    
    disp(['LR voxel size = ',num2str(e1_FOV_mm_SVS),'mm'])
    disp(['AP voxel size = ',num2str(e2_FOV_mm_SVS),'mm'])
    disp(['CC voxel size = ',num2str(e3_FOV_mm_SVS),'mm'])
    
elseif orientation_SVS==2 %coronal
    %     LPS_SVS_centre(1)=-SV_info.rdb_hdr.roilocy;%L
    %     LPS_SVS_centre(2)=-SV_info.rdb_hdr.roilocz;%P
    %     LPS_SVS_centre(3)=SV_info.rdb_hdr.roilocx;%S
    
    e1_FOV_mm_SVS=SV_info.image.user8;%LR
    e2_FOV_mm_SVS=SV_info.image.user10;%CC
    e3_FOV_mm_SVS=SV_info.image.user9;%AP
    
    disp(['LR voxel size = ',num2str(e1_FOV_mm_SVS),'mm'])
    disp(['AP voxel size = ',num2str(e3_FOV_mm_SVS),'mm'])
    disp(['CC voxel size = ',num2str(e2_FOV_mm_SVS),'mm'])
elseif orientation_SVS==1 %sagital
    %     LPS_SVS_centre(1)=-SV_info.rdb_hdr.roilocz;%L
    %     LPS_SVS_centre(2)=-SV_info.rdb_hdr.roilocy;%P
    %     LPS_SVS_centre(3)=SV_info.rdb_hdr.roilocx;%S
    
    e1_FOV_mm_SVS=SV_info.image.user9;%AP
    e2_FOV_mm_SVS=SV_info.image.user10;%CC
    e3_FOV_mm_SVS=SV_info.image.user8;%LR
    
    disp(['LR voxel size = ',num2str(e3_FOV_mm_SVS),'mm'])
    disp(['AP voxel size = ',num2str(e1_FOV_mm_SVS),'mm'])
    disp(['CC voxel size = ',num2str(e2_FOV_mm_SVS),'mm'])
    %     warndlg('MRS volume is not axial. Please double check correct orientation');
end


a=SV_info.rdb_hdr.scan_date;
o=regexp(a,'/','split');
mm=str2double(o{1});
dd=str2double(o{2});
yy=str2double(o{3});
yy=yy+1900;
disp(['Scan date: ',num2str(dd),'/',num2str(mm),'/',num2str(yy)])
disp(['Scan time: ',SV_info.rdb_hdr.scan_time])
if length(num2str(SV_info.rdb_hdr.run_int))==4
    disp(['P0',num2str(SV_info.rdb_hdr.run_int)])
else
    disp(['P',num2str(SV_info.rdb_hdr.run_int)])
end
disp('######################')
end

function [orientation]=check_orient(e3)
[~,orientation]=max(abs(e3));
% orientation==1 --> sagital
% orientation==2 --> coronal
% orientation==3 --> axial
if orientation==1
    disp('sagital volume')
elseif orientation==2
    disp('coronal volume')
elseif orientation==3
    disp('axial volume')
end
end
