% APB: Script to processes the MEGA data
%% Display parameters

function specreg_proc_MEGA(scanner, fmt)

% Process MEGA data    
%
% ARGS :
% scanner = scanner used for acquisition, can be either 
%   Notts-Philips (default)
%   Notts-GE
%   Newcastle
%
% fmt = format of the data, can be 
%   sdat
%   datalist (default)
%   sinlabraw
%
% If no options are given, the default values are used
% RETURNS:
% fids = data after inversed Fourier transform 
%
% EXAMPLE: 
% >> FIDs = mrs_ifft(spectra); 
% >> plot(FIDs)
%
% AUTHOR : Chen Chen
% PLACE  : Sir Peter Mansfield Magnetic Resonance Centre (SPMMRC)

if(nargin<2)
    fmt = 'datalist';
end

if(nargin<1)
    scanner = 'Notts-Philips';
end

%% Processing parameters

spec_lb = 3;
spec_filt = 0.12;
delta0 = 4.7; % is the 'default'

do_DAS          = 1; % causes baseline shifts when large water peak - not recommended for GE
fMRS_block_size = 0; % split into 8 blocks

%% Select the RAW files
flag_Newcastle              = 0;

switch scanner
    case 'Notts-Philips'
        if(strcmp(fmt,'sdat'))
            [data,water,info,filename]  = read_sdat();
        elseif(strcmp(fmt, 'sinlabraw'))
            [data, water, info, filename] = read_sinlabraw();
        elseif(strcmp(fmt, 'datalist'))
            [data,water,info,filename] = read_datalist();
        else
            disp('Unexpected format type')
        end
    case 'Notts-GE'
        [data, water, info, filename] = read_GE();
    case 'Newcastle'
        flag_Newcastle = 1;
        [data,water,info,filename] = read_sdat();
    otherwise
        warning('Error occured when choosing scanner, choose either \n Notts-Philips, Notts-GE, Newcastle')
end

info_struct.fmt         = fmt;
info_struct.scanner     = scanner;
info_struct.filename    = filename;
info_struct.transmit_f  = info.transmit_frequency;
info_struct.TE          = info.TE;
info_struct.TR          = info.TR;
info_struct.B0          = round(abs(info.transmit_frequency/42.57e6));

%% Data parameters

[N, NT]             =   size(data);
NPair               =   NT/2;

[pathname,name,~] = fileparts(filename);
C = strsplit(name,'_');

%% make folders
name=['SpecReg/Data/' C{1}]; % name of your output files [need to adjust accordingly!]
disp(name)

filepath = 'SpecReg/MEGA';
% make a new SpecReg directory

mkdir('SpecReg');

mkdir([filepath '/Figures']);
mkdir([filepath '/Data']);

%% Do an ECC correction
% do it on the data
data = mrs_ifft(data); % convert to FID
data = spa_eddyCor2(mrs_ifft(water(:,2)),data); % Use one of the water refs
data = mrs_fft(data); % convert back


%% Split ON/OFF pairs
if(flag_Newcastle)
    % sometimes requires a -ve flag for some reason?
    data_on         =   mrs_ifft(data(:,2:2:end));
    data_off        =   mrs_ifft(data(:,1:2:end));
else
    data_on         =   mrs_ifft(data(:,1:2:end));
    data_off        =   mrs_ifft(data(:,2:2:end));
end

metab.info      = info;
% Define a 'template' for spectral registration
template        =  mean(data_off,2); % Mean of OFF data

% Option to manually phase the template 1st - shouldn't be needed
[template ,pc] = mrs_manualzeroPHC(mrs_fft(template),[1 N],zeros(1,N));

% Recorrect all spectra with manual phase if applied
data_off    =   mrs_rephase(mrs_fft(data_off),pc);
data_on     =   mrs_rephase(mrs_fft(data_on), pc);

% Convert all spec to FID time domain for specReg
data_off    =   mrs_ifft(data_off);
data_on     =   mrs_ifft(data_on);
template    =   mrs_ifft(template);

%[metab.off_global, f_vec_off, p_vec_off,f_align_OFF]    = spec_reg_fn(data_off, template, metab.info, [0 5], 1, delta0);
% Do first spectral registration to template
[metab.off_global, f_vec_off, p_vec_off,f_align_OFF]    = spec_reg_fn(data_off, template, metab.info, [0 5], 1, delta0);
[metab.on_global, f_vec_on, p_vec_on, f_align_ON]       = spec_reg_fn(data_on, template, metab.info, [0 5], 1, delta0);

figure(f_align_OFF);
print([filepath '/Figures/align_OFF.pdf'], '-dpdf', '-fillpage');
figure(f_align_ON);
print([filepath '/Figures/align_ON.pdf'], '-dpdf', '-fillpage');

% Get correction amount - store in f_vec
f_vec(1:2:size(f_vec_off,2)*2) = f_vec_off;
f_vec(2:2:size(f_vec_off,2)*2) = f_vec_on;
p_vec(1:2:size(f_vec_off,2)*2) = p_vec_off;
p_vec(2:2:size(f_vec_off,2)*2) = p_vec_on;

figure('Name', 'Corrections Applied');
subplot(2,1,1)
plot(f_vec);xlabel('measurement'); ylabel('freq (Hz)'); axis tight; hold on
plot(1:length(f_vec),zeros(1,length(f_vec)), 'r--');legend('F0', 'Template F0');
subplot(2,1,2)
plot(rad2deg(p_vec));xlabel('measurement'); ylabel('phase (deg)');
print([filepath '/Figures/corrections.pdf'], '-dpdf', '-fillpage');

% Have a look at the phase and frequency corrections between ON and OFF

f_vec_diff = abs(f_vec(1:2:end) - f_vec(2:2:end));
p_vec_diff = p_vec(1:2:end) - p_vec(2:2:end);

zfvec = (f_vec_diff - mean(f_vec_diff)) ./std(f_vec_diff);
zpvec = (p_vec_diff - mean(p_vec_diff)) ./std(p_vec_diff);

f_rej = find(zfvec>3);
p_rej = find(abs(zpvec)>3);
if(~isempty(f_rej))
    disp(['Removing shots: ' num2str(f_rej) ' due to extreme (3 sigma) frequency shift of ' num2str(f_vec_diff(f_rej))])
end
if(~isempty(p_rej))
    disp(['Removing shots: ' num2str(p_rej) ' due to extreme (3 sigma) phase shift'])
end

if(sum(f_vec_diff>5))
    disp('WARNING: Frequency shifted > 5 Hz between shots - potential for large drift (FWHM(95%) of 7 T editing pulse = 30Hz)')
end
disp(['Max F shift = ' num2str(max(f_vec_diff))])
disp(['Mean F shift = ' num2str(mean(f_vec_diff))])
disp(['Median F shift = ' num2str(median(f_vec_diff))])
% calculate ppm vector from information
ppm_vec = ppmscale(metab.info.BW, data_on, -metab.info.transmit_frequency/10^6, delta0); % 4.7 because in previous script but need to get exact value

% Get alignment parameters for creatine peak in off-spectrum

%[metab, SNR] = spec_align_to_Cr(metab,ppm_vec)

%% Shift NAA frequency to actual ppm units

[metab, SNR] = spec_reference_NAA(metab,ppm_vec);

info_struct.SNR = SNR;

%% now align pairs
for n= 1:size(data_off,2)
    [data_off_algn(:,n), f_vec(n), p_vec(n)] = spec_reg_fn(metab.off_global(:,n), metab.on_global(:,n), metab.info, [2.8 3.15], 1, delta0, 0, 0.2);
end

%data_diff_pair = metab.on_global-data_off_algn;

%for n= 1:size(data_diff_pair,2)
%    data_diff_algn(:,n) = spec_reg_fn(data_diff_pair(:,n), mean(data_diff_pair,2), metab.info, [0 5], 1, delta0, 0.2);
%end

metab.off_global = data_off_algn;



%% Automatic rejection of datasets

% plot before rejections
f=figure('Name', 'AutoReject Check');plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(mean(metab.on_global,2)-mean(metab.off_global,2),spec_lb,info.BW,spec_filt))));
set(gca, 'XDir', 'reverse');xlabel('ppm');

% choose Cho peak to assess outlier spectra (3.1-3.3)ppm
P               = (ppm_vec > 3.1 & ppm_vec < 3.3);
metab.concat(:,1:2:size(f_vec_off,2)*2)    = metab.off_global;
metab.concat(:,2:2:size(f_vec_off,2)*2)    = metab.on_global;

% Get real-valued spectrum around Cho peak only = R
R = mrs_fft(metab.concat);
R = real(R(P,:)); % real cho peak

% For all measurements - measure the mean square error MSE from mean
% -could try to zerofill to help here

% try zerofill

for i = 1:size(R,2)
    MSE(i) = norm(mean(R,2)-R(:,i));
end
% compute z-score
zMSE = (MSE - mean(MSE)) ./std(MSE);

% threshold around central value
Q1 = prctile(MSE, 25);
Q3 = prctile(MSE, 75);
% Compute the IQR
IQR = Q3 - Q1;
% Define outlier thresholds
lower_bound = Q1 - 1.5 * IQR;
upper_bound = Q3 + 1.5 * IQR;

figure('Name', 'Automatic Outlier Rejection');
subplot(1,2,1)
rej_idx = find(abs(zMSE)>3);
plot(zMSE); title('MSE/mean(MSE) (Cho peak)');
ylabel('z-stat, units of standard deviation');
xlabel('measurement');
hold on; scatter(rej_idx, zMSE(rej_idx),200, 'xr', 'LineWidth', 2);
legend('MSE', 'Rejected- 3 sigma', 'Location', 'southwest' );
subplot(1,2,2)
plot(ppm_vec, real(mean(mrs_fft(metab.concat),2)), 'k', 'LineWidth', 1.2);
hold on
labelRej{1} = 'Mean';


if(~isempty(rej_idx))
disp('****')
    for r = 1:length(rej_idx)
        disp(['Removing shot ' num2str(rej_idx(r)) ' z-stat: ' num2str(zMSE(rej_idx(r)))])
        plot(ppm_vec, real(mrs_fft(broaden_filter_FID_sw(metab.concat(:,r),0, info.BW, 0.10))));
        labelRej{r+1} = num2str(rej_idx(r));
    end
else
    disp('No shots removed')
end

set(gca, 'XDir', 'reverse');
xlim([2.4 4.2]);
legend(labelRej, 'Location', 'northwest' );
xlabel('ppm');
title('Rejected Spectra');
print([filepath '/Figures/rejections.pdf'], '-dpdf', '-fillpage');

% remove these and then average (=MEAN) all the spectra ON and OFF
OFF_rej = [];
ON_rej  = [];
for r=1:length(rej_idx)
    if(mod(rej_idx(r),2)) % this is odd = OFF
        OFF_rej = [OFF_rej 1+(rej_idx(r)-1)/2];
    else
        ON_rej = [ON_rej rej_idx(r)/2];
    end
end

metab.off_rejected = metab.off_global;
metab.on_rejected = metab.on_global;

% Throw away bad shots
% sometimes there are duplicates - need to ensure throwing right ones away

rej_array = unique([OFF_rej, ON_rej, f_rej, p_rej]); % find unique indices to throw

metab.off_rejected(:, rej_array) = [];
%metab.off_rejected(:, rej_array) = [];
metab.on_rejected(:, rej_array) = []; % make sure to now throw out pairs
%metab.on_rejected(:, rej_array) = [];


%% Water referencing
% use the water collected in same file 
% but for newcastle data use NWS

if(flag_Newcastle)
    disp('choose the NWS file');
    [water,~,info,filename] = read_sdat(); % get NWS
    % Do ECC of water reference
    waterf = mrs_ifft(water); % convert to FID
    waterf = spa_eddyCor2(waterf(:,1),waterf); % Eddy current corretion
else
    % Do ECC of water reference
    waterf = mrs_ifft(water); % convert to FID
    waterf = spa_eddyCor2(waterf(:,2),waterf); % Eddy current corretion
end

% need to correct water references
[water_aligned, f_vec, p_vec, f_align]    = spec_reg_fn(waterf, waterf(:,1), metab.info, [3 6.5], 0, delta0);
waterfid = mean(water_aligned,2); % average water acquisitions
txfrq = info.transmit_frequency;

[LW] = meas_LW_water(waterf, metab,ppm_vec);

info_struct.LW = LW;

%% fMRS split here
if(fMRS_block_size > 0) % if splitting into Blocks
    mkdir([filepath '/Data/Block/']);
    block_path = [filepath '/Data/Block/' num2str(fMRS_block_size)];
    mkdir(block_path);
    
    NBlocks = NPair / fMRS_block_size;

    for nB = 1:NBlocks
        tmp_block = (1+(fMRS_block_size)*(nB-1):(fMRS_block_size)*nB);
        %metab.off_global_blocks(:,:,nB) = metab.off_global(:,tmp_block);
        %metab.on_global_blocks(:,:,nB) = metab.on_global(:,tmp_block);
    
    %% do SR in blocks
    template        =  mean(metab.off_global(:,tmp_block),2);
    [metab.off_global_blocks(:,:,nB), f_vec_off, p_vec_off,f_align_OFF]    = spec_reg_fn(metab.off_global(:,tmp_block), template, metab.info, [0 5], 1, delta0);
    [metab.on_global_blocks(:,:,nB), f_vec_on, p_vec_on, f_align_ON]       = spec_reg_fn(metab.on_global(:,tmp_block), template, metab.info, [0 5], 1, delta0);
    end
    
    % remove the outliers from the averaged spectra

    metab.off_rejected_blocks = metab.off_global_blocks;
    metab.on_rejected_blocks = metab.on_global_blocks;

    [row, col] = ind2sub([fMRS_block_size, NBlocks], OFF_rej);
    metab.off_rejected_blocks(:, row, col) = NaN; % make NAN

    [row, col] = ind2sub([fMRS_block_size, NBlocks], ON_rej);
    metab.on_rejected_blocks(:, row, col) = NaN; % make NAN

    s1 = squeeze(nanmean(metab.off_rejected_blocks,2));
    s2 = squeeze(nanmean(metab.on_rejected_blocks,2));

    s2_cor = run_DAS(s1, s2, metab, delta0, 0, filepath); % DAS alignment - OFF for fMRS
    
    metab   = save_param(metab, s1,s2_cor);
    save_diff_mat(s1, s2_cor, block_path);
    
    write_sdatspar(metab,filename,[block_path '/' C{1}]);
    write_spafiles(metab, waterfid, txfrq, info, delta0, [block_path '/' C{1}]);
end

%% Main Alignment

s1=squeeze(mean(metab.off_rejected,2));
s2=squeeze(mean(metab.on_rejected,2));


s2_cor = run_DAS(s1, s2, metab, delta0, do_DAS, filepath);

% diff_res = metab.on_rejected - metab.off_rejected;
% 
% for n= 1:size(diff_res,2)
%    diff_res_c(:,n) = spec_reg_fn(diff_res(:,n), mean(diff_res,2), metab.info, [2.75 3.25], 1, delta0, 0, 0.1);
% end

metab   = save_param(metab, s1 ,s2_cor);
save_diff_mat(s1, s2_cor, [filepath '/Data']);

write_sdatspar(metab,filename,[filepath '/Data/' C{1}]);
write_spafiles(metab, waterfid, txfrq, info, delta0, [filepath '/Data/' C{1}]);


% save information
info_struct.name = C{1};
JSONFILE_name= sprintf([filepath '/info.json']); 
fid=fopen(JSONFILE_name,'w');
encodedJSON = jsonencode(info_struct);
fprintf(fid, encodedJSON);
fclose(fid);

%% Plots

figure('Name', 'Final ON/OFF before subtraction');
plot(ppm_vec, real(mrs_fft(s2_cor)));
hold on
plot(ppm_vec, real(mrs_fft(s1)));
legend('ON', 'OFF');
set(gca, 'XDir', 'reverse');
xlim([0.5 4.2]);
xlabel('ppm');
print([filepath '/Figures/final_OnOff.pdf'], '-dpdf', '-fillpage');

figure('Name', 'SR comparison');hold on
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(s2-s1,spec_lb,info.BW,spec_filt))));
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(s2_cor-s1,spec_lb,info.BW,spec_filt))));
yL = get_ylim_editing(real(mrs_fft(broaden_filter_FID_sw(s2_cor-s1,spec_lb,info.BW,spec_filt))), ppm_vec);
set(gca, 'XDir', 'reverse');
xlim([0 5]);
ylim(yL);
xlabel('ppm');
legend('SR','SR+DAS');
print([filepath '/Figures/SRDAScomparison.pdf'], '-dpdf', '-fillpage');

% plot after rejections
figure(f);hold on;plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(s2-s1,spec_lb,info.BW,spec_filt))));
set(gca, 'XDir', 'reverse');title('Diff Spec');xlabel('ppm');
xlim([2.2 4.2]);legend('Before Rejections', 'After Rejections');
print([filepath '/Figures/rejections_influence.pdf'], '-dpdf');

figure('Name', 'No Corrections');
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(sum(data_on,2)-sum(data_off,2),spec_lb,info.BW,spec_filt))));
yL = get_ylim_editing(real(mrs_fft(broaden_filter_FID_sw(sum(data_on,2)-sum(data_off,2),spec_lb,info.BW,spec_filt))), ppm_vec);
set(gca, 'XDir', 'reverse');
xlim([0 5]);
ylim(yL);
xlabel('ppm');
legend('No Corr');
print([filepath '/Figures/spec_no_corr.pdf'], '-dpdf', '-fillpage');

figure('Name', 'Diff Spec');hold on
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw((s2_cor-s1)*exp(1i*0*pi/180),spec_lb,info.BW,spec_filt))));
%hold on
%plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(mean(diff_res_c,2),spec_lb,info.BW,spec_filt))));
yL=get_ylim_editing(real(mrs_fft(broaden_filter_FID_sw((s2_cor-s1)*exp(1i*0*pi/180),spec_lb,info.BW,spec_filt))), ppm_vec);
set(gca, 'XDir', 'reverse');
xlim([0 5]);
ylim(yL);
xlabel('ppm');
print([filepath '/Figures/difference_spectrum.pdf'], '-dpdf', '-fillpage');


figure('Name', 'Diff Spec _ scaled');hold on
pbaspect([1.0000    0.2046    0.2046])
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw((s2_cor-s1)*exp(1i*0*pi/180),spec_lb,info.BW,spec_filt))));
%hold on
%plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(mean(diff_res,2),spec_lb,info.BW,spec_filt))));
yL=get_ylim_editing(real(mrs_fft(broaden_filter_FID_sw((s2_cor-s1)*exp(1i*0*pi/180),spec_lb,info.BW,spec_filt))), ppm_vec);
set(gca, 'XDir', 'reverse');
xlim([0 5]);
ylim(yL);
xlabel('ppm');
print([filepath '/Figures/difference_spectrum_scaled.pdf'], '-dpdf', '-fillpage');

end

%% Functions
function s2_cor = run_DAS(s1, s2, metab, delta0, do_DAS, filepath)

if(do_DAS)
    for i = 1:size(s2,2)
             [s2_cor, f_vec_local, ph_vec_local,f_align_DAS] = spec_reg_fn(s2(:,i), s1(:,i), metab.info, [3.1 3.3],1,delta0); % previously 3.15-3.35
    end
    figure(f_align_DAS);
    print([filepath '/Figures/align_DAS.pdf'], '-dpdf', '-fillpage');
else
    s2_cor = s2;
end

end


function metab_out = save_param(metab_in, s1,s2_cor)
    metab_out = metab_in;
    metab_out.final.on      = s2_cor;
    metab_out.final.off     = s1;
    metab_out.final.diff    = s2_cor-s1;
    metab_out.final.sum     = s2_cor+s1;
end

function save_diff_mat(s1, s2_cor, filepath)
    diff_FID    = s2_cor-s1;
    sum_FID     = s2_cor+s1;
    save([filepath '/diff_spec.mat'], 'diff_FID');
    save([filepath '/sum_spec.mat'], 'sum_FID');
end

function write_sdatspar(metab, filename, filepath)

    outpath = './';

    mrs_writeSDAT([outpath,filepath,'_on.SDAT'], mrs_fft(metab.final.on));
    mrs_writeSDAT([outpath,filepath,'_off.SDAT'], mrs_fft(metab.final.off));
    mrs_writeSDAT([outpath,filepath,'_diff.SDAT'], mrs_fft(metab.final.diff));
    mrs_writeSDAT([outpath,filepath,'_sum.SDAT'], mrs_fft(metab.final.sum));

    copyfile(filename,[outpath,filepath,'_on.SPAR']);
    copyfile(filename,[outpath,filepath,'_off.SPAR']);
    copyfile(filename,[outpath,filepath,'_diff.SPAR']);
    copyfile(filename,[outpath,filepath,'_sum.SPAR']);

end

function write_spafiles(metab, waterfid, txfrq, info, delta0, filepath)

    outpath = './';
    
    write_spa([outpath,filepath,'_diff.spa'], (metab.final.diff), waterfid, txfrq, info.BW, delta0);
    write_spa([outpath,filepath,'_off.spa'], (metab.final.off), waterfid, txfrq, info.BW, delta0);
    write_spa([outpath,filepath,'_on.spa'], (metab.final.on), waterfid, txfrq, info.BW, delta0);
    write_spa([outpath,filepath,'_sum.spa'], (metab.final.sum), waterfid, txfrq, info.BW, delta0);

end

