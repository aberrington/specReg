function [metab, SNR] = spec_reference_NAA(metab,ppm_vec)


%3.15 3.4
%spec_peak = 3.185
% between 0 and 4 ppm - look for max peak height - this will be NAA (at 
tmp=find((ppm_vec>0 & ppm_vec<4));
rng = [tmp(1) tmp(end)];
[A I] = mrs_findPeak(mrs_fft(sum(metab.off_global,2)), rng);
rng_ppm_NAA = 0.8; % ppm
ppm_per_point = ppm_vec(2)-ppm_vec(1);
shft = round(rng_ppm_NAA/ppm_per_point/2);
[A_peak, I_peak, peak_fitted, pars_fitted] = mrs_fitPeak( mrs_fft(sum(metab.off_global,2)), [(I-abs(shft)),(I + abs(shft))], 0, 'l'); % fit peak
I_peak_ppm = ppm_vec(I_peak);

spec_shift_ppm = (2.01-I_peak_ppm);
disp(['Shifting spectra by ' num2str(spec_shift_ppm) ' ppm']);
spec_shift = round((2.01-I_peak_ppm)/ppm_per_point);
test = circshift((sum(metab.off_global,2)),-spec_shift);

for i = 1:size(metab.off_global,2)
    metab.off_global(:,i)   = mrs_ifft(circshift(mrs_fft(metab.off_global(:,i)),spec_shift));
    metab.on_global(:,i)    = mrs_ifft(circshift(mrs_fft(metab.on_global(:,i)),spec_shift));
end

in.ppm = ppm_vec;
in.specs = conj(mrs_fft(squeeze(nanmean(metab.off_global,2))));

[SNR, signal,noisesd]=op_getSNR(in,1.5,2.2,-10,-5); % get SNR of OFF peak of NAA


end

