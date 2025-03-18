function [fids_aligned] = spec_reg_fn_apply(fids, f_vec, ph_vec, info, freq_range, isCut, delta0, cutDur)
% This script is going to do spectral registration algorithm of all fids to
% a template. Can limit this over a certain range (ppm).
% based on Near et al. (2015)
% INPUT
% fids: this is the time domain FID data
% template: this is the template to align to
% need to restrict frequency firstly fft, then cut, then ifft
if(nargin<8)
    cutDur = 0.2; % 200ms as default
end
if(nargin<7)
    delta0 = 4.7;
end
if(nargin<6)
    isCut = 1;
end
if(nargin<5)
    freq_range = [3.15 3.35];
end
%freq_range = [-2 4.2];
dt      = 1/info.BW;
t2 = dt*(0:(size(fids,1)-1))';

%template_zf = mrs_zerofill(template, 0);
fids_zf  = mrs_zerofill(fids, 0);
if(isCut)
    cut_spec = cutDur/dt;
else
    cut_spec = size(template,1);
end
%template_filt = broaden_filter_FID_sw(template_zf, 0,info.BW, Inf);
%template_filt(0.2/dt:end,:)=0; % null all above 200ms
%template_filt = template_filt(1:cut_spec); % null all above 200ms
fids_filt = broaden_filter_FID_sw(fids_zf, 0,info.BW, Inf);
%fids_filt(0.2/dt:end,:)  = 0; % null all above 200ms
fids_filt = fids_filt(1:cut_spec,:); % null all above 200ms

%ppm_vec = ppmscale(info.BW, fids_filt, -info.transmit_frequency/10^6, 4.7); % 4.7 because in previous script but need to get exact value
% may have to remove some samples at end of FID for 'smoothing', also
% linebroadening and zeropad?
%template_f  = mrs_fft(template_filt);
ppm_vec = ppmscale(info.BW, fids_filt, -info.transmit_frequency/10^6, delta0); % 4.7 because in previous script but need to get exact value
P           = (ppm_vec > freq_range(1) & ppm_vec < freq_range(2));

%template_f_cut = template_f(P);
%T = mrs_ifft(template_f_cut); % this is the template over the freq range
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

for j = 1:size(fids,2)
    fids_aligned(:,j) = fids(:,j) .* exp(1i*2*pi*(f_vec(j).*t2 + ph_vec(j)));
end

% get time-domain representation over freq_range

% 
% figure;
% hold on
% plot(real(mrs_fft(T))); % template
% plot(real(mrs_fft(F(:,30)))); % orig
% plot(real(mrs_fft(F(:,30).*exp(1i*2*pi*(f_vec(30).*t + ph_vec(30)))))); % reg
% %params = -params;

% 
end
