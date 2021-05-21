% APB_read_GE_Notts
function [data, water, info, filename] = read_GE()

zf_factor=1;

%% read data
if nargin==0
    [filename,filepath]=uigetfile('*.7', 'Select GE Pfile','/home/acrad/Rdrive/Academic-Radiology');
end

header=read_MR_headers([filepath,filename],'all');
%% MRS acquisition info
disp('-----------')
disp([filename,':'])
get_MRS_info(header);

SV_info.num_points=header.rdb_hdr.frame_size; 
SV_info.bw_Hz=header.rdb_hdr.user0;
% SV_info.nex=1;
SV_info.nex=8;
% % SV_info.nex=header.image.nex;
SV_info.chop = header.rdb_hdr.user41;                            % chop due to no_add data
if (SV_info.chop)
    SV_info.nsa  = header.rdb_hdr.user4;                        % with no_add option
else    
    SV_info.nsa  = header.rdb_hdr.user4/header.rdb_hdr.navs;    % data frames
end
SV_info.coils=(header.rdb_hdr.dab(2)-header.rdb_hdr.dab(1))+1;
SV_info.transmit_frequency = -header.rdb_hdr.ps_mps_freq/10;
ppm=(linspace(SV_info.bw_Hz/2000,-SV_info.bw_Hz/2000,zf_factor*SV_info.num_points)/(128*10^3))*10^6+4.7;

% don't need below - hard coded! 
% fid = fopen([filepath,filename],'r', 'ieee-le');
% status = fseek(fid, 1468, 'bof');
% pfile_header_size = fread(fid,1,'integer*4');

pfile_header_size = header.total_length; % APB this is already in header

fid = fopen([filepath,filename],'r', 'ieee-le');
status = fseek(fid, pfile_header_size, 'bof');
raw_data=fread(fid,'int32');


num_spectra=length(raw_data)/(SV_info.num_points*2);
all_spectra=zeros(zf_factor*SV_info.num_points,num_spectra);

fid_processed=zeros(zf_factor*SV_info.num_points,1);

x=linspace(0,(zf_factor*SV_info.num_points-1)/SV_info.bw_Hz,zf_factor*SV_info.num_points)';
lb_Hz=0; %optional exp line broadening in Hz
apo1=exp(-x*lb_Hz);


for t=1:num_spectra
    offset=SV_info.num_points*2*(t-1);
    fid_re=raw_data(offset+1:2:offset+SV_info.num_points*2);
    fid_im=raw_data(offset+2:2:offset+SV_info.num_points*2);
    
    fid=fid_re+1i*fid_im;
    
    fid_processed(1:SV_info.num_points,1)=fid;    
    fid_processed(:,1)=fid_processed(:,1).*apo1;   
    
    all_spectra(:,t)=fftshift(fft(fid_processed(:,1)))./1e6;%some random scaling to keep numbers small
end

% clear empty frames
all_spectra(:,1:header.rdb_hdr.nframes+1:end)=[];

clear raw_data
clear fid_processed
clear fid

%% sort data into water (WS) and metabo (nonWS)
% B is unedited and A is edited (might be the other way around depending on
% scanner settings)
counter=0;
for coil=1:2:(SV_info.coils*2)  
    counter=counter+1;
    water_A(:,:,counter)=all_spectra(:,1+(header.rdb_hdr.nframes*(coil-1)):SV_info.nex+(header.rdb_hdr.nframes*(coil-1)));
    water_B(:,:,counter)=all_spectra(:,1+(header.rdb_hdr.nframes*(coil)):SV_info.nex+(header.rdb_hdr.nframes*(coil)));    
    metabo_A(:,:,counter)=all_spectra(:,SV_info.nex+1+(header.rdb_hdr.nframes*(coil-1)):SV_info.nsa+SV_info.nex+(header.rdb_hdr.nframes*(coil-1)));
    metabo_B(:,:,counter)=all_spectra(:,SV_info.nex+1+(header.rdb_hdr.nframes*(coil)):SV_info.nsa+SV_info.nex+(header.rdb_hdr.nframes*(coil)));       
end

% although the nonWS is stored sparately before the WS A and B data, there
% are no editing pulses applied to water_A, hence putting water_A and
% water_B together
water=[water_A water_B];

clear all_spectra
%% Klose ECC correction=0 order phasing --> only on water averages thus not very accurate but looks like method of choice
% set1 and set2 are used because there is a 180degree phase shift in the data   
water_set1_coil=squeeze(mean(water(:,1:2:end,:),2));
water_set2_coil=squeeze(mean(water(:,2:2:end,:),2));
metabo_A_set1=metabo_A(:,1:2:end,:);
metabo_A_set2=metabo_A(:,2:2:end,:);

metabo_B_set1=metabo_B(:,1:2:end,:);
metabo_B_set2=metabo_B(:,2:2:end,:);

%do ECC
for t=1:(SV_info.nsa/2)
    for coil=1:SV_info.coils   
         metabo_phase_A_set1(:,t,coil) = fftshift(fft(ECCKlose(ifft(ifftshift(metabo_A_set1(:,t,coil))),ifft(ifftshift(water_set1_coil(:,coil))))));
         metabo_phase_A_set2(:,t,coil) = fftshift(fft(ECCKlose(ifft(ifftshift(metabo_A_set2(:,t,coil))),ifft(ifftshift(water_set2_coil(:,coil))))));
         
         metabo_phase_B_set1(:,t,coil) = fftshift(fft(ECCKlose(ifft(ifftshift(metabo_B_set1(:,t,coil))),ifft(ifftshift(water_set1_coil(:,coil))))));
         metabo_phase_B_set2(:,t,coil) = fftshift(fft(ECCKlose(ifft(ifftshift(metabo_B_set2(:,t,coil))),ifft(ifftshift(water_set2_coil(:,coil))))));         
    end
end


%ECC eliminates phase shift, so can join the two sets
metabo_phase_A=[metabo_phase_A_set1 metabo_phase_A_set2];
clear metabo_phase_A_set1
clear metabo_phase_A_set2

metabo_phase_B=[metabo_phase_B_set1 metabo_phase_B_set2];
clear metabo_phase_B_set1
clear metabo_phase_B_set2




%% get coil weightings from unsupressed water and avg water signal for qunatification
for t=1:16
    if mod(t,2)
        [phase_ref(:,t),water_fp(:,t,:)]=first_point_phase(squeeze(water(:,t,:)));
    else
        [phase_ref(:,t),water_fp(:,t,:)]=first_point_phase(squeeze(-water(:,t,:)));
    end
end
water_mean=squeeze(mean(water_fp,2));
coil_w=squeeze(max(real(water_mean)));
water_mean = water_mean * coil_w';
for nsa_w = 1:size(water_fp,2)
    
    water_shots(:,nsa_w)=squeeze(water_fp(:,nsa_w,:))*coil_w';

end



%% coil combination
for nsa=1:SV_info.nsa
    dummy=squeeze(metabo_phase_A(:,nsa,:));
    all_A_nsa(:,nsa)=dummy*coil_w';  
    dummy=squeeze(metabo_phase_B(:,nsa,:));
    all_B_nsa(:,nsa)=dummy*coil_w';  
    
    SNR(nsa)=get_SNR(ppm,all_A_nsa(:,nsa));
end
[~,m]=max(SNR);

all_A=squeeze(mean(all_A_nsa,2));
all_B=squeeze(mean(all_B_nsa,2));
%write_RAW_LCModel([filepath,filename, 'GABAprocess_ON.RAW'],(ifft(ifftshift(all_A))),filename,header)
%write_RAW_LCModel([filepath,filename, 'GABAprocess_OFF.RAW'],(ifft(ifftshift(all_B))),filename,header)
%write_RAW_LCModel([filepath,filename, 'GABAprocess_ON-OFF.RAW'],(ifft(ifftshift(all_A-all_B))),filename,header)
%write_RAW_LCModel([filepath,filename, 'GABAprocess_water_ref.RAW'],(ifft(ifftshift(water_mean))),filename,header)

% [water_out,spec_out]=water_removal_felix(all_A-all_B,4096);

figure;
subplot(2,3,1);plot(ppm,real(all_A_nsa));set(gca,'xdir','reverse');xlim([0 4]);title('A')
subplot(2,3,2);plot(ppm,real(all_B_nsa));set(gca,'xdir','reverse');xlim([0 4]);title('B')
subplot(2,3,4);plot(ppm,real(all_A));set(gca,'xdir','reverse');xlim([0 4]);title('A')
subplot(2,3,5);plot(ppm,real(all_B));set(gca,'xdir','reverse');xlim([0 4]);title('B')
subplot(2,3,3);plot(ppm,real(all_A-all_B));set(gca,'xdir','reverse');xlim([0 4]);title('A-B')

lb_on_off=ifft(ifftshift(all_A-all_B));
x=linspace(0,(zf_factor*SV_info.num_points-1)/SV_info.bw_Hz,zf_factor*SV_info.num_points)';
lb_Hz=5; %optional exp line broadening in Hz
apo1=exp(-x*lb_Hz);
lb_on_off=lb_on_off.*apo1;
lb_on_off=fftshift(fft(lb_on_off));
subplot(2,3,6);plot(ppm,real(lb_on_off));set(gca,'xdir','reverse');xlim([0 4]);title('A-B lb=5Hz')

saveas(gcf,[filepath,filename,'_GABAprocess_specREGPRE.png'],'png')

water   = water_shots;
data = zeros(SV_info.num_points,nsa*2);
data(:,1:2:end)         = all_A_nsa;
data(:,2:2:end)         = all_B_nsa;
info.BW     = SV_info.bw_Hz;
info.transmit_frequency = SV_info.transmit_frequency;
