function [LW] = meas_LW_water(waterf, metab,ppm_vec)

tmp=find((ppm_vec>3 & ppm_vec<6));
rng = [tmp(1) tmp(end)];
[A I] = mrs_findPeak(mrs_fft(sum(waterf,2)), rng);
rng_ppm_NAA = 4; % ppm
ppm_per_point = ppm_vec(2)-ppm_vec(1);
shft = round(rng_ppm_NAA/ppm_per_point/2);
[A_peak, I_peak, peak_fitted, pars_fitted] = mrs_fitPeak( mrs_fft(sum(waterf,2)), [(I-abs(shft)),(I + abs(shft))], 0, 'l'); % fit peak
LW = pars_fitted(3)*metab.info.BW/length(ppm_vec);

disp(['The estimated water linewidth is: ' num2str(LW) 'Hz'])
end

