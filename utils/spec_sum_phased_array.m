function [outfidw,outfidu,weights,outfid_wref]= spec_sum_phased_array(infid,inref,mrprot,quiet,plot_result)
% sum spectra acquired using phased array coil with phase correction

% first, zero-pad and gaussian apodize
% expfilt = exp(-(1:2048).*(1:2048)/(384^2));
% inref = [inref ; 0*inref];
% infid = [infid ; 0*infid];
% for x=1:8
%     inref(:,x) = inref(:,x) .* expfilt';
%     infid(:,x) = infid(:,x) .* expfilt';
% end

[dx dz] = size(inref);

%correct eddy currents and 0-order phase for all spectra
% phi_w = angle(inref);
% phi_w_ref = phi_w(1,:);
% phi_m = angle(infid);
% phi_m_ref = phi_m(1,:);
% 
% phi_w = phi_w - ones(dx,1)*phi_w_ref;
% %phi_w_ref = phi_w(1,:);
% phi_m = phi_m - ones(dx,1)*phi_m_ref - phi_w;
% 
% 
% abs_fid = abs(infid);
% abs_ref = abs(inref);
% 
% pcorref = complex(abs_ref .* cos(phi_w), abs_ref .* sin(phi_w));
% pcorfid = complex(abs_fid .* cos(phi_m), abs_fid .* sin(phi_m));

pcorref = inref;
pcorfid = infid;
for x=1:dz
    pcorref(:,x) = phase_adjust(inref(:,x),angle(inref(1,x)));
    pcorfid(:,x) = phase_adjust(infid(:,x),angle(inref(1,x)));
end

% next, compute weights by fitting a high-order polynomial and taking the scale
% factor
xdata = (1:dx)';
warning off MATLAB:polyfit:RepeatedPointsOrRescale
pcoef = zeros(dz,10);  % 9th order (9+1)
for x=1:dz
    pcoef(x,:) = polyfit(xdata,abs(double(pcorref(:,x))),9); % 9th order
end
weights = pcoef(:,end) / sum(pcoef(:,end));
[null,sortix] = sort(weights,'descend');
if plot_result
 plot8(real(fftshift(fft(double(inref(:,sortix))),1)),'uncor ref =',1:8);
 plot8(real(fftshift(fft(double(pcorref(:,sortix))),1)),'cor ref =',1:8);
 plot8(real(fftshift(fft(double(infid(:,sortix))),1)),'uncor met =',1:8);
 plot8(real(fftshift(fft(double(pcorfid(:,sortix))),1)),'cor met =',1:8);
end

% compute weighted and unweighted sums
unweightsum = 0 * pcorfid;
weightsum = 0 * pcorfid;
unweightsum_wref=0*pcorref;
weightsum_wref=0*pcorref;

for x=1:dz
    unweightsum(:,x:dz) = unweightsum(:,x:dz) + repmat(1/dz * pcorfid(:,sortix(x)),1,dz-x+1);
    unweightsum_wref(:,x:dz) = unweightsum_wref(:,x:dz) + repmat(1/dz * pcorref(:,sortix(x)),1,dz-x+1);
    weightsum(:,x:dz) = weightsum(:,x:dz) + repmat(weights(sortix(x)) * pcorfid(:,sortix(x)),1,dz-x+1);
    weightsum_wref(:,x:dz) = weightsum_wref(:,x:dz) + repmat(weights(sortix(x)) * pcorref(:,sortix(x)),1,dz-x+1);
end

% compute signal based on region from 1.5-2.5 ppm (naa, hopefully)
sigstart = ppmpx(mrprot, dx, 2.5);
sigfinish = ppmpx(mrprot, dx, 1.5);

% tmppcorfid = real(fftshift(fft(double(pcorfid(:,sortix))),1))
% tmppcorfid = tmppcorfid(sigstart:sigfinish,:);
% plot8(tmppcorfid,'cor naa only =',1:8);

% compute noise sd based on region from 8-11 ppm
sdstart = ppmpx(mrprot, dx, 11);
if (sdstart < 50), sdstart = 50; end  % don't look at the very edge
sdfinish = ppmpx(mrprot, dx, 8);

% singleft = real(fftshift(fft(pcorfid),1));
% singlemax = max(singleft(sigstart:sigfinish,:));
% singlenoise = std(singleft(sdstart:sdfinish,:));
% 
% unweightft = real(fftshift(fft(unweightsum),1));
% unweightmax = max(unweightft(sigstart:sigfinish,:));
% unweightnoise = std(unweightft(sdstart:sdfinish,:));
% 
% weightft = real(fftshift(fft(weightsum),1));
% weightmax = max(weightft(sigstart:sigfinish,:));
% weightnoise = std(weightft(sdstart:sdfinish,:));

if ~quiet
    fprintf('\n');
    for x=1:dz
        fprintf('single coil (#%2d): max = %10.2f ; sd = %10.2f ; snr = %6.2f ; weight = %4.1f%%\n',sortix(x),singlemax(sortix(x)),singlenoise(sortix(x)),singlemax(sortix(x))/singlenoise(sortix(x)),100*weights(sortix(x)));
    end
    fprintf('\n');
    for x=1:dz
        fprintf('unweight sum 1-%2d: max = %10.2f ; sd = %10.2f ; snr = %6.2f\n',x,unweightmax(x),unweightnoise(x),unweightmax(x)/unweightnoise(x));
    end
    fprintf('\n');
    for x=1:dz
        fprintf('weighted sum 1-%2d: max = %10.2f ; sd = %10.2f ; snr = %6.2f\n',x,weightmax(x),weightnoise(x),weightmax(x)/weightnoise(x));
    end
    fprintf('\n');

    fprintf('* Weights are computed based on unsuppressed water signal.\n')
    fprintf('* Signal estimates are based on assumed NAA peak (maximum between 1.5-2.5 ppm).\n')
    fprintf('* Noise estimate is SD of the region between 8-11 ppm.\n')
    fprintf('* Weights and SNR will not necessarily correlate, but should be close.\n')

    fprintf('\n');
end
if plot_result
figure; plot(singleft(:,sortix(1))); axis([1 dx (-.5 * max(singleft(:,sortix(1)))) 1.1*max(singleft(:,sortix(1)))])
title('Best single-coil spectrum')
figure; plot(unweightft(:,end)); axis([1 dx (-.5 * max(unweightft(:,end))) 1.1*max(unweightft(:,end))])
title('Unweighted summed spectrum')
figure; plot(weightft(:,end)); axis([1 dx (-.5 * max(weightft(:,end))) 1.1*max(weightft(:,end))])
title('Weighted summed spectrum')
end
outfidw = weightsum(:,end);
outfidu = unweightsum(:,end);
outfid_wref=weightsum_wref(:,end);


