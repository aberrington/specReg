% APB: Script to processes the MEGA-sLASER data for 7 T data using Raw
% format
%% Display parameters

spec_lb = 3;
spec_filt = 99;
delta0 = 4.7; % is the 'default'

%% Processing parameters

do_DAS          = 1; % causes baseline shifts when large water peak - not recommended for GE
fMRS_block_size = 41; % split into 8 blocks

%% Select the RAW files
% This is for 7 T Raw files
%[data,water,info,filename] = read_sinlabraw();
[data, water, info, filename] = read_datalist();

% This is for 3 T Philips
%[data,water,info,filename] = read_sdat();
flag_Newcastle = 0; % Newcastle = 1

% This is for GE (work in progress)
%[data, water, info, filename] = read_GE();

% so far not processed the NWS water references for Newc
% some Newcastle data require a flip of everything ...
%% Data parameters

[N, NT]             =   size(data);
NPair               =   NT/2;

[pathname,name,~] = fileparts(filename);
C = strsplit(name,'_');

%% make folders
name=['SpecReg/Data/' C{1}]; % name of your output files [need to adjust accordingly!]
disp(name)

% make a new SpecReg directory

mkdir('SpecReg');
mkdir('SpecReg/Figures');
mkdir('SpecReg/Data');


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

% Do first spectral registration to template
[metab.off_global, f_vec_off, p_vec_off,f_align_OFF]    = spec_reg_fn(data_off, template, metab.info, [0 5], 1, delta0);
[metab.on_global, f_vec_on, p_vec_on, f_align_ON]       = spec_reg_fn(data_on, template, metab.info, [0 5], 1, delta0);

figure(f_align_OFF);
print('SpecReg/Figures/align_OFF.pdf', '-dpdf', '-fillpage');
figure(f_align_ON);
print('SpecReg/Figures/align_ON.pdf', '-dpdf', '-fillpage');

% Get correction amount - store in f_vec
f_vec(1:2:size(f_vec_off,2)*2) = f_vec_off;
f_vec(2:2:size(f_vec_off,2)*2) = f_vec_on;
p_vec(1:2:size(f_vec_off,2)*2) = p_vec_off;
p_vec(1:2:size(f_vec_off,2)*2) = p_vec_on;

figure('Name', 'Corrections Applied');
subplot(2,1,1)
plot(f_vec);xlabel('measurement'); ylabel('freq (Hz)'); axis tight; hold on
plot(1:length(f_vec),zeros(1,length(f_vec)), 'r--');legend('F0', 'Template F0');
subplot(2,1,2)
plot(rad2deg(p_vec));xlabel('measurement'); ylabel('phase (deg)');
print('SpecReg/Figures/corrections.pdf', '-dpdf', '-fillpage');
% calculate ppm vector from information
ppm_vec = ppmscale(metab.info.BW, data_on, -metab.info.transmit_frequency/10^6, delta0); % 4.7 because in previous script but need to get exact value

%% Shift NAA frequency to actual ppm units

[metab, SNR] = spec_reference_NAA(metab,ppm_vec);

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
for i = 1:size(R,2)
    MSE(i) = norm(mean(R,2)-R(:,i));
end
% compute z-score
zMSE = (MSE - mean(MSE)) ./std(MSE);
 
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
        disp(['Removing shot ' num2str(rej_idx(r)) ' z-stat: ' num2str(zMSE(rej_idx))])
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
print('SpecReg/Figures/rejections.pdf', '-dpdf', '-fillpage');

% remove these and then average (=MEAN) all the spectra ON and OFF
OFF_rej = [];
ON_rej  = [];
for r=1:length(rej_idx)
    if(mod(rej_idx(r),2)) % this is odd = OFF
        OFF_rej = 1+(rej_idx(r)-1)/2;
    else
        ON_rej = rej_idx(r)/2;
    end
end

metab.off_rejected = metab.off_global;
metab.on_rejected = metab.on_global;

% Throw away bad shots
metab.off_rejected(:, OFF_rej) = [];
metab.on_rejected(:, ON_rej) = [];


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

%% fMRS split here
if(fMRS_block_size > 0) % if splitting into Blocks
    mkdir('SpecReg/Data/Block/');
    block_path = ['SpecReg/Data/Block/' num2str(fMRS_block_size)];
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

    s2_cor = run_DAS(s1, s2, metab, delta0, 0); % DAS alignment - OFF for fMRS
    
    metab   = save_param(metab, s1,s2_cor);
    save_diff_mat(s1, s2_cor, block_path);
    
    write_sdatspar(metab,filename,[block_path '/' C{1}]);
    write_spafiles(metab, waterfid, txfrq, info, delta0, [block_path '/' C{1}]);
end

%% Main Alignment

s1=squeeze(mean(metab.off_rejected,2));
s2=squeeze(mean(metab.on_rejected,2));

s2_cor = run_DAS(s1, s2, metab, delta0, do_DAS);

metab   = save_param(metab, s1 ,s2_cor);
save_diff_mat(s1, s2_cor, 'SpecReg/Data');

write_sdatspar(metab,filename,['SpecReg/Data/' C{1}]);
write_spafiles(metab, waterfid, txfrq, info, delta0, ['SpecReg/Data/' C{1}]);

%% Plots

figure('Name', 'Final ON/OFF before subtraction');
plot(ppm_vec, real(mrs_fft(s2_cor)));
hold on
plot(ppm_vec, real(mrs_fft(s1)));
legend('ON', 'OFF');
set(gca, 'XDir', 'reverse');
xlim([0.5 4.2]);
xlabel('ppm');
print('SpecReg/Figures/final_OnOff.pdf', '-dpdf', '-fillpage');

figure('Name', 'SR comparison');hold on
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(s2-s1,spec_lb,info.BW,spec_filt))));
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(s2_cor-s1,spec_lb,info.BW,spec_filt))));
yL = get_ylim_editing(real(mrs_fft(broaden_filter_FID_sw(s2_cor-s1,spec_lb,info.BW,spec_filt))), ppm_vec);
set(gca, 'XDir', 'reverse');
xlim([0 5]);
ylim(yL);
xlabel('ppm');
legend('SR','SR+DAS');
print('SpecReg/Figures/SRDAScomparison.pdf', '-dpdf', '-fillpage');

% plot after rejections
figure(f);hold on;plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(s2-s1,spec_lb,info.BW,spec_filt))));
set(gca, 'XDir', 'reverse');title('Diff Spec');xlabel('ppm');
xlim([2.2 4.2]);legend('Before Rejections', 'After Rejections');
print('SpecReg/Figures/rejections_influence.pdf', '-dpdf');

figure('Name', 'No Corrections');
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(sum(data_on,2)-sum(data_off,2),spec_lb,info.BW,spec_filt))));
yL = get_ylim_editing(real(mrs_fft(broaden_filter_FID_sw(sum(data_on,2)-sum(data_off,2),spec_lb,info.BW,spec_filt))), ppm_vec);
set(gca, 'XDir', 'reverse');
xlim([0 5]);
ylim(yL);
xlabel('ppm');
legend('No Corr');
print('SpecReg/Figures/spec_no_corr.pdf', '-dpdf', '-fillpage');

figure('Name', 'Diff Spec');hold on
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw((s2_cor-s1)*exp(1i*0*pi/180),spec_lb,info.BW,spec_filt))));
yL=get_ylim_editing(real(mrs_fft(broaden_filter_FID_sw((s2_cor-s1)*exp(1i*0*pi/180),spec_lb,info.BW,spec_filt))), ppm_vec);
set(gca, 'XDir', 'reverse');
xlim([0 5]);
ylim(yL);
xlabel('ppm');
print('SpecReg/Figures/difference_spectrum.pdf', '-dpdf', '-fillpage');


%% Functions
function s2_cor = run_DAS(s1, s2, metab, delta0, do_DAS)

if(do_DAS)
    for i = 1:size(s2,2)
             [s2_cor, f_vec_local, ph_vec_local,f_align_DAS] = spec_reg_fn(s2(:,i), s1(:,i), metab.info, [3.1 3.3],1,delta0); % previously 3.15-3.35
    end
    figure(f_align_DAS);
    print('SpecReg/Figures/align_DAS.pdf', '-dpdf', '-fillpage');
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

