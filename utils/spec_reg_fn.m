function [fids_aligned, f_vec, ph_vec, fig_align] = spec_reg_fn(fids, template, info, freq_range, isCut, delta0, cutDur)
% This script is going to do spectral registration algorithm of all fids to
% a template. Can limit this over a certain range (ppm).
% based on Near et al. (2015)
% INPUT
% fids: this is the time domain FID data
% template: this is the template to align to
% need to restrict frequency firstly fft, then cut, then ifft
if(nargin<7)
    cutDur = 0.2; % 200ms as default
end
if(nargin<6)
    delta0 = 4.7;
end
if(nargin<5)
    isCut = 1;
end
if(nargin<4)
    freq_range = [3.15 3.35];
end
%freq_range = [-2 4.2];
dt      = 1/info.BW;
t2 = dt*(0:(size(fids,1)-1))';

template_zf = mrs_zerofill(template, 0);
fids_zf  = mrs_zerofill(fids, 0);
if(isCut)
    cut_spec = cutDur/dt;
else
    cut_spec = size(template,1);
end
template_filt = broaden_filter_FID_sw(template_zf, 0,info.BW, Inf);
%template_filt(0.2/dt:end,:)=0; % null all above 200ms
template_filt = template_filt(1:cut_spec); % null all above 200ms
fids_filt = broaden_filter_FID_sw(fids_zf, 0,info.BW, Inf);
%fids_filt(0.2/dt:end,:)  = 0; % null all above 200ms
fids_filt = fids_filt(1:cut_spec,:); % null all above 200ms

%ppm_vec = ppmscale(info.BW, fids_filt, -info.transmit_frequency/10^6, 4.7); % 4.7 because in previous script but need to get exact value
% may have to remove some samples at end of FID for 'smoothing', also
% linebroadening and zeropad?


% get time-domain representation over freq_range
template_f  = mrs_fft(template_filt);
ppm_vec = ppmscale(info.BW, fids_filt, -info.transmit_frequency/10^6, delta0); % 4.7 because in previous script but need to get exact value
P           = (ppm_vec > freq_range(1) & ppm_vec < freq_range(2));

template_f_cut = template_f(P);
T = mrs_ifft(template_f_cut); % this is the template over the freq range
% need new time vector from this
t = dt*(0:(sum(P==1)-1))'; % this is the time vector

fid_f = mrs_fft(fids_filt);
fid_f_cut = fid_f(P,:);
F = mrs_ifft(fid_f_cut);

params0(1) = 0; % initial 0 freq shift
params0(2) = 0; % initial 0 phase shift

%lb = [zeros(8,1); -1 * Inf(8,1)]; % set the lower bound to the solution - amplitudes should be < 1 and > 0
%ub = [ones(8,1); Inf(8,1)];
%[params, fval] = fmincon(@(params)(costMinTargetAmp(params(:),S)),params0,[],[],[],[],[], [],@limit_channel_RF_norm,options);

cutSamples = length(find(P(:)==1));

%BWperSample     = info.BW/size(template,1);
%factor          = cutSamples/(cutSamples*BWperSample);
time_conv_factor = (cutSamples/cut_spec);
for j = 1:size(F,2)
    [params,fval] = fminsearch(@(params)(costSpectralReg(params(:),T, F(:,j), t)),params0);
    f_vec(j) = params(1);
    ph_vec(j) = params(2);
    fids_aligned(:,j) = fids(:,j) .* exp(1i*2*pi*(params(1).*t2*time_conv_factor + params(2)));
end

% 
% figure;
% hold on
% plot(real(mrs_fft(T))); % template
% plot(real(mrs_fft(F(:,30)))); % orig
% plot(real(mrs_fft(F(:,30).*exp(1i*2*pi*(f_vec(30).*t + ph_vec(30)))))); % reg
% %params = -params;

% 
if(size(F,2)==1)
    fig_align = figure('Name', 'Spec Reg Alignment');
    hold on
    plot(real(mrs_fft(T))); % template
    plot(real(mrs_fft(F))); % orig
    plot(real(mrs_fft(F(:,1:end).*exp(1i*2*pi*(f_vec(1:end).*t + ph_vec(1:end)))))); % reg
    legend('temp', 'orig', 'SR');
else
    fig_align = figure('Name', 'Spec Reg Alignment');
    subplot(2,1,1)
    plot(real(mrs_fft(F))); % orig
    title('Before');
    subplot(2,1,2)
    plot(real(mrs_fft(F(:,1:end).*exp(1i*2*pi*(f_vec(1:end).*t + ph_vec(1:end)))))); % reg
    title('After');
end
% % 
% figure;
% hold on
% plot(real(mrs_fft(template))); % template
% plot(real(mrs_fft(fids(:,1)))); % orig
% plot(real(mrs_fft(fids(:,1).*exp(1i*2*pi*(params(1)*t2*time_conv_factor + params(2)))))); % reg
% legend('template', 'orig', 'SR');
% %params = -params;
f_vec = f_vec *time_conv_factor;

end

function  cost = costSpectralReg(params,template,fid,t)
% similar to MRCode tools option 5

    f   = params(1);
    phi = params(2);
    
    G = fid .* exp(1i*2*pi*(f*t + phi));
    
    cost = norm(template - G);
end