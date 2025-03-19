function y = mrs_lorentzFun_new( x, y0, x0, fwhm, A  )
% MRS_LORENTZFUN defines the lorentzian function.  
% 
% y = mrs_lorentzFun( x, y0, x0, fwhm, A  )
%
% ARGS :
% x = input vectors 
% y0 = baseline amplitude
% x0 = location of the peak 
% fwhm = full width at half maximum
% A = height of the peak 
% 
% RETURNS:
% y = output vectors
%
%
% AUTHOR : Chen Chen - modified by A Berrington
% PLACE  : Sir Peter Mansfield Magnetic Resonance Centre (SPMMRC)
%
     gamma = fwhm/2;

     y = y0 + (A/pi) * (gamma / ((x-x0).^2 + gamma^2);


     %m=(x-x0)*2./fwhm;
     %y= y0 + A.*(1./(1+m.^2));
     
end

