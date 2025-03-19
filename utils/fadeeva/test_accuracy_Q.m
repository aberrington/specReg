% Set 34 digits, so IEEE quadruple precision, in the ADVANPIX 
% Multiprecision Computing Toolbox (advanpix.com)
mp.Digits(34)
% Create a 2001X801 array containing the 1,602,801 complex numbers
% z = 10^p exp(i theta), with p = -6(0.0006)6, theta = 0(pi/1600)pi/2
p = -6:0.0006:6;
theta = (0:0.1125:90)*pi/180;
z = zeros(length(p),length(theta));
for i = 1: length(p)
    for j = 1: length(theta)
        z(i,j) = 10^p(i)*exp(1i*theta(j));
    end
end
% Calculate accurate values for w(z), by the call wTrap_Q(z,N)
% which uses our new method with N = 20 evaluated in quadruple 
% precision via the ADVANPIX toolbox.
truew = double(wTrap_Q(z,20)); wabs = abs(truew);
% Compute maximum absolute and relative errors over these z values for ...
% ... our new method with N = 11
errors = abs(truew-wTrap(z,11)); 
a1 = max(max(errors)); r1 = max(max(errors./wabs));
% ... Weideman's approximation with N = 40
errors = abs(truew-cef(z,40)); 
a2 = max(max(errors)); r2 = max(max(errors./wabs));
% ... Abrarov et al's approximation
errors = abs(truew-fadsamp(z)); 
a3 = max(max(errors)); r3 = max(max(errors./wabs));
% ... Zaghloul and Ali's approximation
errors = abs(truew-Faddeyeva_v2(z,13)); 
a4 = max(max(errors)); r4 = max(max(errors./wabs));
% Display these errors
format SHORTG 
Absolute = [a1;a2;a3;a4]
Relative = [r1;r2;r3;r4]