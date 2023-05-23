% APB: Script to processes the single shot data
%% Display parameters

spec_lb = 3;
spec_filt = 0.12;
delta0 = 4.7; % is the 'default'


%% Select the RAW files
% This is for 7 T Raw files
%[data,water,info,filename] = read_sinlabraw();
[data, water, info, filename] = read_datalist();

% This is for 3 T Philips
%[data,water,info,filename] = read_sdat();

% This is for GE (work in progress)  
%[data, water, info, filename] = read_GE();

% so far not processed the NWS water references for Newc
% some Newcastle data require a flip of everything ...
%% Data parameters

[N, NT]             =   size(data);

[pathname,name,~] = fileparts(filename);
C = strsplit(name,'_');

%% make folders
name=['SpecReg/SLASER/Data/' C{1}]; % name of your output files [need to adjust accordingly!]
disp(name)

% make a new SpecReg directory
filepath = 'SpecReg/SLASER';

mkdir('SpecReg');

mkdir([filepath '/Figures']);
mkdir([filepath '/Data']);


%% Do an ECC correction
% do it on the data
data = mrs_ifft(data); % convert to FID
data = spa_eddyCor2(mrs_ifft(water(:,1)),data); % Use one of the water refs
data = mrs_fft(data); % convert back


%% Split ON/OFF pairs

    % sometimes requires a -ve flag for some reason?
dataC            =   mrs_ifft(data);
metab.info       =   info;

% Define a 'template' for spectral registration
template        =  mean(dataC,2); % Mean of OFF data

% Option to manually phase the template 1st - shouldn't be needed
[template ,pc] = mrs_manualzeroPHC(mrs_fft(template),[1 N],zeros(1,N));

% Recorrect all spectra with manual phase if applied
dataC    =   mrs_rephase(mrs_fft(dataC),pc);

% Convert all spec to FID time domain for specReg
dataC        =   mrs_ifft(dataC);
template    =   mrs_ifft(template);

% Do first spectral registration to template
%metab.info.transmit_frequency = -metab.info.transmit_frequency
[metab.global, f_vec, p_vec, f_align]    = spec_reg_fn(dataC, template, metab.info, [-Inf Inf], 1, 4.65, 0.05);
% use first 100 ms of FID
%fidCor = spa_xcorr(dataC, info.BW);
%[fidCor,phzmethod] = spa_pcorr(fidCor,1, info.BW, 4.65, info.transmit_frequency/1e6);

figure
plot(real(mean(mrs_fft(metab.global),2)))
hold on
plot(real(mean(mrs_fft(dataC),2)))
legend('SR', 'orig')

figure(f_align);
print([filepath '/Figures/aligned.pdf'], '-dpdf', '-fillpage');

figure('Name', 'Corrections Applied');
subplot(2,1,1)
plot(f_vec);xlabel('measurement'); ylabel('freq (Hz)'); axis tight; hold on
plot(1:length(f_vec),zeros(1,length(f_vec)), 'r--');legend('F0', 'Template F0');
subplot(2,1,2)
plot(rad2deg(p_vec));xlabel('measurement'); ylabel('phase (deg)');
print([filepath '/Figures/corrections.pdf'], '-dpdf', '-fillpage');
% calculate ppm vector from information
ppm_vec = ppmscale(metab.info.BW, data, -metab.info.transmit_frequency/10^6, delta0); % 4.7 because in previous script but need to get exact value

%% Shift NAA frequency to actual ppm units

[metab, SNR] = spec_reference_NAA(metab,ppm_vec,0,0);

%% Automatic rejection of datasets

% plot before rejections
f=figure('Name', 'AutoReject Check');plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw(mean(metab.global,2),spec_lb,info.BW,spec_filt))));
set(gca, 'XDir', 'reverse');xlabel('ppm');

% choose Cho peak to assess outlier spectra (3.1-3.3)ppm
P               = (ppm_vec > 3.1 & ppm_vec < 3.3);

% Get real-valued spectrum around Cho peak only = R
R = mrs_fft(metab.global);
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
plot(ppm_vec, real(mean(mrs_fft(metab.global),2)), 'k', 'LineWidth', 1.2);
hold on
labelRej{1} = 'Mean';


if(~isempty(rej_idx))
disp('****')
    for r = 1:length(rej_idx)
        disp(['Removing shot ' num2str(rej_idx(r)) ' z-stat: ' num2str(zMSE(rej_idx(r)))])
        plot(ppm_vec, real(mrs_fft(broaden_filter_FID_sw(metab.global(:,r),0, info.BW, 0.10))));
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


metab.global_rejected = metab.global;

% Throw away bad shots
metab.global_rejected(:, rej_idx) = [];


%% Water referencing
% use the water collected in same file 
% but for newcastle data use NWS

    % Do ECC of water reference
waterf = mrs_ifft(water); % convert to FID
waterf = spa_eddyCor2(waterf(:,1),waterf); % Eddy current corretion


% need to correct water references
[water_aligned, f_vec, p_vec, f_align]    = spec_reg_fn(waterf, waterf(:,1), metab.info, [3 6.5], 1, delta0);
waterfid = mean(water_aligned,2); % average water acquisitions
txfrq = info.transmit_frequency;

[LW] = meas_LW_water(waterf, metab,ppm_vec);


%% Main Alignment

s=squeeze(mean(metab.global_rejected,2));

metab   = save_param(metab, s);
save_data(s, [filepath '/Data']);

write_sdatspar(metab,filename,[filepath '/Data/' C{1}]);
write_spafiles(metab, waterfid, txfrq, info, delta0, [filepath '/Data/' C{1}]);

%% Plots


figure('Name', 'Spec');hold on
plot(ppm_vec,real(mrs_fft(broaden_filter_FID_sw((s)*exp(1i*0*pi/180),spec_lb,info.BW,spec_filt))));
yL=get_ylim_editing(real(mrs_fft(broaden_filter_FID_sw((s)*exp(1i*0*pi/180),spec_lb,info.BW,spec_filt))), ppm_vec);
set(gca, 'XDir', 'reverse');
xlim([0 5]);
%ylim(yL);
xlabel('ppm');
print([filepath '/Figures/spectrum.pdf'], '-dpdf', '-fillpage');



%% Functions

function metab_out = save_param(metab_in, s)
    metab_out = metab_in;
    metab_out.final      = s;
end

function save_data(s1, filepath)
    FID    = s1;
    save([filepath '/spec.mat'], 'FID');
    
end

function write_sdatspar(metab, filename, filepath)

    outpath = './';

    mrs_writeSDAT([outpath,filepath,'.SDAT'], mrs_fft(metab.final));
    copyfile(filename,[outpath,filepath,'.SPAR']);

end

function write_spafiles(metab, waterfid, txfrq, info, delta0, filepath)

    outpath = './';
    write_spa([outpath,filepath,'.spa'], metab.final, waterfid, txfrq, info.BW, delta0);

end

