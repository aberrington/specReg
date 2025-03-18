function [metab] = spec_align_to_Cr(metab,ppm_vec,isedit,isshift)

if(nargin <4)
    isshift=1;
end

if(nargin<3)
    isedit=1;
end

if(~isedit)
    metab.off_global = metab.global;
end
%3.15 3.4
%spec_peak = 3.185
% between 0 and 4 ppm - look for max peak height - this will be Cr
tmp=find((ppm_vec>2.72 & ppm_vec<3.12));
rng = [tmp(1) tmp(end)];
[A I] = mrs_findPeak(mrs_fft(sum(metab.off_global,2)), rng);
rng_ppm_NAA = 0.2; % ppm
ppm_per_point = ppm_vec(2)-ppm_vec(1);
shft = round(rng_ppm_NAA/ppm_per_point/2);

for i= 1:size(metab.off_global,2)
    [A_peak, I_peak, peak_fitted, pars_fitted] = mrs_fitPeak_phase( mrs_fft(metab.off_global(:,i)), [(I-abs(shft)),(I + abs(shft))], 0, 'lp'); % fit peak
    I_peak_ppm = ppm_vec(I_peak);
    spec_shift_ppm = (3.027-I_peak_ppm);
    disp(['Shifting spectra by ' num2str(spec_shift_ppm) ' ppm']);
    spec_shift = round((3.027-I_peak_ppm)/ppm_per_point);
    test = circshift((sum(metab.off_global,2)),-spec_shift);
    
    if(isshift)
        metab.off_global(:,i)   = mrs_ifft(circshift(mrs_fft(metab.off_global(:,i)),spec_shift));
        metab.off_global(:,i)   = metab.off_global(:,i)*exp(-1i*pars_fitted(5)); % apply phase too        
        metab.on_global(:,i)    = mrs_ifft(circshift(mrs_fft(metab.on_global(:,i)),spec_shift));
        metab.on_global(:,i)    = metab.on_global(:,i)*exp(-1i*pars_fitted(5)); % apply phase too
     else
         disp('Not shifting the spectra')
     end
end


% 
% in.ppm = ppm_vec';
% in.specs = conj(mrs_fft(squeeze(mean(metab.off_global,2, 'omitnan'))));
% 
% [SNR, signal,noisesd]=op_getSNR_specReg(in,1.5,2.2); % get SNR of OFF peak of NAA

if(~isedit)
    metab.global = metab.off_global;
end

end

