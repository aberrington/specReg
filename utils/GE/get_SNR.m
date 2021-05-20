function [SNR]=get_SNR(ppm,spec_in)
% SNR=max(signal:4ppm-0ppm)/std(-6ppm:end)


%Determine limits
z=abs(ppm-4.0);
lb=find(min(z)==z);
z=abs(ppm-0.0);
ub=find(min(z)==z);

z=abs(ppm+6.0);
noise=find(min(z)==z);


SNR=max(real(spec_in(lb:ub)))/std(real(spec_in(noise:end)));



end