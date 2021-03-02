% APB: Script to processes the MEGA-sLASER data for 7 T data using Raw
% format

%% Select the RAW files
[filename1,~,~] = uigetfile('.raw','Select .raw file');
[filename2,pathname,~] = uigetfile('.SPAR','Select .SPAR file');
cd(pathname);

% make a new SpecReg directory
mkdir('SpecReg');
mkdir('SpecReg/Figures');
mkdir('SpecReg/Data');

[~,name,~] = fileparts(filename2);
C = strsplit(name,'_');

name=['SpecReg/Data/' C{3}]; % name of your output files [need to adjust accordingly!]
disp(name)

[data, water]       =   mrs_readRAW(filename1); % reads in the sin/lab/raw data and gets the metab/water specs
info                =   mrs_readSPAR(filename2);
[N, NT]             =   size(data);
NPair               =   NT/2;

%% Do an ECC correction
data = mrs_ifft(data); % convert to FID
data = spa_eddyCor2(mrs_ifft(water(:,2)),data); % Use one of the water refs
data = mrs_fft(data); % convert back

%% Split ON/OFF pairs
data_on         =   mrs_ifft(data(:,1:2:end));
data_off        =   mrs_ifft(data(:,2:2:end));

metab.info      = info;
% Define a 'template' for spectral registration
template        =  mean(data_off,2); % Mean of OFF data

% Option to manually phase the template 1st - shouldn't be needed
[template ,pc] = mrs_manualzeroPHC(mrs_fft(template),[1 4096],zeros(1,4096));

% Recorrect all spectra with manual phase if applied
data_off    =   mrs_rephase(mrs_fft(data_off),pc);
data_on     =   mrs_rephase(mrs_fft(data_on), pc);

% Convert all spec to FID time domain for specReg
data_off    =   mrs_ifft(data_off);
data_on     =   mrs_ifft(data_on);
template    =   mrs_ifft(template);

% Do first spectral registration to template
[metab.off_global, f_vec_off, p_vec_off,f_align_OFF]    = spec_reg_fn(data_off, template, metab.info, [0 5]);
[metab.on_global, f_vec_on, p_vec_on, f_align_ON]       = spec_reg_fn(data_on, template, metab.info, [0 5]);

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
ppm_vec = ppmscale(metab.info.BW, data_on, -metab.info.transmit_frequency/10^6, 4.7); % 4.7 because in previous script but need to get exact value


%% Automatic rejection of datasets

% plot before rejections
f=figure('Name', 'AutoReject Check');plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(mean(metab.on_global,2)-mean(metab.off_global,2),0,info.BW,0.10))));
set(gca, 'XDir', 'reverse');xlabel('ppm');

% choose Cho peak to assess outlier spectra (3.15-3.35)ppm
P               = (ppm_vec > 3.15 & ppm_vec < 3.35);
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
        OFF_rej = (rej_idx(r)-1)/2;
    else
        ON_rej = rej_idx(r)/2;
    end
end

metab.off_rejected = metab.off_global;
metab.on_rejected = metab.on_global;

% Throw away bad shots
metab.off_rejected(:, OFF_rej) = [];
metab.on_rejected(:, ON_rej) = [];

s1=squeeze(mean(metab.off_rejected,2));
s2=squeeze(mean(metab.on_rejected,2));

% plot after rejections
figure(f);hold on;plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(s2-s1,0,info.BW,0.10))));
set(gca, 'XDir', 'reverse');title('Diff Spec');xlabel('ppm');
xlim([2.2 4.2]);legend('Before Rejections', 'After Rejections');
print('SpecReg/Figures/rejections_influence.pdf', '-dpdf');
%% DAS - Final SR of ON and OFF before subtraction

%Do spectralReg over Cho peak only ... 3.15 to 3.35
for i = 1:size(s2,2)
     [s2_cor, f_vec_local, ph_vec_local,f_align_DAS] = spec_reg_fn(s2, s1, metab.info, [3.15 3.35],1);
end
figure(f_align_DAS);
print('SpecReg/Figures/align_DAS.pdf', '-dpdf', '-fillpage');

figure('Name', 'Final ON/OFF before subtraction');
plot(ppm_vec, real(mrs_fft(s2_cor)));
hold on
plot(ppm_vec, real(mrs_fft(s1)));
legend('ON', 'OFF');
set(gca, 'XDir', 'reverse');
xlim([0.5 4.2]);
xlabel('ppm');
print('SpecReg/Figures/final_OnOff.pdf', '-dpdf', '-fillpage');


figure('Name', 'No Corrections');
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(sum(data_on,2)-sum(data_off,2),3,info.BW,0.10))));
yL = get_ylim_editing(real(mrs_fft(broaden_filter_FID_sw(sum(data_on,2)-sum(data_off,2),3,info.BW,0.10))), ppm_vec)
set(gca, 'XDir', 'reverse');
xlim([0 5]);
ylim(yL);
xlabel('ppm');
legend('No Corr');
print('SpecReg/Figures/spec_no_corr.pdf', '-dpdf', '-fillpage');

figure('Name', 'SR comparison');hold on
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(s2-s1,3,info.BW,0.10))));
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(s2_cor-s1,3,info.BW,0.10))));
yL = get_ylim_editing(real(mrs_fft(broaden_filter_FID_sw(s2_cor-s1,3,info.BW,0.10))), ppm_vec);
set(gca, 'XDir', 'reverse');
xlim([0 5]);
ylim(yL);
xlabel('ppm');
legend('SR','SR+DAS');
print('SpecReg/Figures/SRDAScomparison.pdf', '-dpdf', '-fillpage');


figure('Name', 'Diff Spec');hold on
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw((s2_cor-s1)*exp(1i*0*pi/180),6,info.BW,0.10))));
yL=get_ylim_editing(real(mrs_fft(broaden_filter_FID_sw((s2_cor-s1)*exp(1i*0*pi/180),3,info.BW,0.10))), ppm_vec);
set(gca, 'XDir', 'reverse');
xlim([0 5]);
ylim(yL);
xlabel('ppm');
print('SpecReg/Figures/difference_spectrum.pdf', '-dpdf', '-fillpage');

metab.final.on      = s2_cor;
metab.final.off     = s1;
metab.final.diff    = s2_cor-s1;
metab.final.sum     = s2_cor+s1;

diff_FID    = s2_cor-s1;
sum_FID     = s2_cor+s1;

save('SpecReg/Data/diff_spec.mat', 'diff_FID');
save('SpecReg/Data/sum_spec.mat', 'sum_FID');
%% Export to SDAT/SPAR as is currently done in the pipeline...

mrs_writeSDAT([pathname, name,'_on.SDAT'], mrs_fft(metab.final.on));
mrs_writeSDAT([pathname,name,'_off.SDAT'], mrs_fft(metab.final.off));
mrs_writeSDAT([pathname,name,'_diff.SDAT'], mrs_fft(metab.final.diff));
mrs_writeSDAT([pathname,name,'_sum.SDAT'], mrs_fft(metab.final.sum));
%mrs_writeSDAT([pathname,name,'_sum.SDAT'], sum_avg);

copyfile(filename2,[pathname,name,'_on.SPAR']);
copyfile(filename2,[pathname,name,'_off.SPAR']);
copyfile(filename2,[pathname,name,'_diff.SPAR']);
copyfile(filename2,[pathname,name,'_sum.SPAR']);
%copyfile(filename2,[pathname,name,'_sum.SPAR']);




