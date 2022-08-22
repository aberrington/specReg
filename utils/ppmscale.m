function scale_ppm=ppmscale(sw,spectra,sfrq1H,H1offset)
if(nargin<4)
    H1offset = 4.65;
end

%H1offset = 0;
fmax = (sw)/2;
f = fmax: -2*fmax/(length(spectra)-1): -fmax;
scale_ppm = f/(sfrq1H) + H1offset;
end
