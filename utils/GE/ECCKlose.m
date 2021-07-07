% Klose's Eddy Current Correction
%signal: complex FID
%watersignal: complex FID of the unsupressed water signal
% This function corrects for eddy currents using Klose's method
function signal = ECCKlose(signal,watersignal)
    phasesignal = angle(signal)-angle(watersignal);
    signal = abs(signal).*cos(phasesignal)+1i.*abs(signal).*sin(phasesignal);

