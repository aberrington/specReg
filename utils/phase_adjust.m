function outfid = phase_adjust(trc, phase)
% function to apply phase angles to complex data, e.g. to use
% variable receiver phases in VNMR 

rtrc = real(trc);
itrc = imag(trc);

cosp = cos(phase); % real_cor
sinp = sin(phase); % imag_cor

rout = rtrc.*cosp + itrc.*sinp;
iout = itrc.*cosp - rtrc.*sinp;

outfid = complex(rout, iout);
