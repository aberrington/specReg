function y = mrs_voigtFun( x, y0, x0, A, sigma, gamma)
% MRS_VOIGTFUN defines the voigt function.  
% 
% y = mrs_voigtFun( x, y0, x0, fwhm, A  )
%
% ARGS :
% x = input vectors 
% y0 = baseline amplitude
% x0 = location of the peak 
% fwhm = full width at half maximum
% A = height of the peak
% gamma = 
% 
% RETURNS:
% y = output vectors
%
%
% AUTHOR : Chen Chen
% PLACE  : Sir Peter Mansfield Magnetic Resonance Centre (SPMMRC)
%
% Copyright (c) 2013, University of Nottingham. All rights reserved.

    z = (x - x0 + sqrt(-1)*gamma) / (sigma * sqrt(2));
    w = wTrapWCP(z,11); % 11 recommended points for eval.

    y = A * real(w) / (sigma * sqrt(2*pi)) + y0;
    %s=fwhm/2.3548; 
    %y=A./s./sqrt(2*pi).*exp(-(x-x0).^2/2/s.^2)+y0;
     
end

