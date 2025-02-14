function [FWHM] = meas_LW_water(waterf, metab,ppm_vec, model)

if(nargin<4)
    model='v';
end
tmp=find((ppm_vec>3 & ppm_vec<6));
rng = [tmp(1) tmp(end)];
[A I] = mrs_findPeak(mrs_fft(sum(waterf,2)), rng);
rng_ppm_NAA = 4; % ppm
ppm_per_point = ppm_vec(2)-ppm_vec(1);
shft = round(rng_ppm_NAA/ppm_per_point/2);
if(strcmp(model, 'v'))
    [A_peak, I_peak, peak_fitted, pars_fitted] = mrs_fitPeak_voigt( mrs_fft(sum(waterf,2)), [(I-abs(shft)),(I + abs(shft))], 0); % fit peak
    FWHM = 0.5346*pars_fitted(5) + sqrt(0.2166*pars_fitted(5)^2 + pars_fitted(3)^2); % using definition of width of voigt profile
elseif(strcmp(model, 'l'))
    [A_peak, I_peak, peak_fitted, pars_fitted] = mrs_fitPeak_lorentz( mrs_fft(sum(waterf,2)), [(I-abs(shft)),(I + abs(shft))], 0); % fit peak
    FWHM = pars_fitted(3);
else
    [A_peak, I_peak, peak_fitted, pars_fitted] = mrs_fitPeak( mrs_fft(sum(waterf,2)), [(I-abs(shft)),(I + abs(shft))], 0, 'g'); % guassian fit peak
    FWHM = pars_fitted(3);
end

FWHM = FWHM*metab.info.BW/length(ppm_vec);

disp(['The estimated water linewidth is: ' num2str(FWHM) 'Hz'])
end

