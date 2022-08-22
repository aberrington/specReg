function [outfidw,outfidu,weights,outfid_wref]= sum_phased_array(infid,inref,quiet,plot_result)
% sum spectra acquired using phased array coil with phase correction

% Code obtained from Uzay Emir
% Edited by Adam Berrington Aug 2017



% Edited by Adam Berrington Aug 2017
[dn, dc]     =   size(inref);

% Use first water reference as the reference for phase and amplitude
% correction

% Then correct others based off this.

pcorref = zeros(dn,dc);
pcorfid = zeros(dn,dc);

for x = 1:dc
    
    pcorref(:,x) = phase_adjust(inref(:,x),angle(inref(1,x)));
    pcorfid(:,x) = phase_adjust(infid(:,x),angle(inref(1,x)));
%     figure
%     plot(1:2048, squeeze(inref(x,1,:)))
%     hold on
%     plot(1:2048, squeeze(pcorref(x,1,:)));
%     pause
end

% next, compute weights by fitting a high-order polynomial and taking the scale
% factor
xdata   = (1:dn)';
warning off MATLAB:polyfit:RepeatedPointsOrRescale
pcoef   = zeros(dc,10);  % 9th order (9+1)

for x = 1:dc
    pcoef(x,:)  = polyfit(xdata,abs(double(squeeze(pcorref(:,x)))),9); % 9th order
end

weights         = pcoef(:,end) / sum(pcoef(:,end));
[~ , sortix]    = sort(weights,'descend');

if plot_result
    figure
    plot(real(fftshift(fft(double(squeeze(inref))),1)))
    figure
    plot(real(fftshift(fft(double(squeeze(pcorref(:,sortix)))),1)))
    figure
    plot(real(fftshift(fft(double(squeeze(infid(:,sortix)))),1)))
    figure
    plot(real(fftshift(fft(double(squeeze(pcorfid(:,sortix)))),1)))
end

% compute weighted and unweighted sums
unweightsum         =   0 * pcorfid;
weightsum           =   0 * pcorfid;
unweightsum_wref    =   0 * pcorref;
weightsum_wref      =   0 * pcorref;

for x=1:dc
    unweightsum(:,x:dc)         = unweightsum(:,x:dc)           + repmat(1/dc * pcorfid(:, sortix(x)),1,dc-x+1);
    unweightsum_wref(:,x:dc)    = unweightsum_wref(:,x:dc)      + repmat(1/dc * pcorref(:, sortix(x)),1,dc-x+1);
    weightsum(:,x:dc)           = weightsum(:,x:dc)             + repmat(weights(sortix(x)) * pcorfid(:, sortix(x)),1,dc-x+1);
    weightsum_wref(:,x:dc)      = weightsum_wref(:,x:dc)        + repmat(weights(sortix(x)) * pcorref(:, sortix(x)),1,dc-x+1);
end

% compute signal based on region from 1.5-2.5 ppm (naa, hopefully)
% sigstart = ppmpx(mrprot, dx, 2.5);
% sigfinish = ppmpx(mrprot, dx, 1.5);

% tmppcorfid = real(fftshift(fft(double(pcorfid(:,sortix))),1))
% tmppcorfid = tmppcorfid(sigstart:sigfinish,:);
% plot8(tmppcorfid,'cor naa only =',1:8);

% % compute noise sd based on region from 8-11 ppm
% sdstart = ppmpx(mrprot, dx, 11);
% if (sdstart < 50), sdstart = 50; end  % don't look at the very edge
% sdfinish = ppmpx(mrprot, dx, 8);

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
% if plot_result
% % figure; plot(singleft(:,sortix(1))); axis([1 dx (-.5 * max(singleft(:,sortix(1)))) 1.1*max(singleft(:,sortix(1)))])
% % title('Best single-coil spectrum')
% % figure; plot(unweightft(:,end)); axis([1 dx (-.5 * max(unweightft(:,end))) 1.1*max(unweightft(:,end))])
% % title('Unweighted summed spectrum')
% % figure; plot(weightft(:,end)); axis([1 dx (-.5 * max(weightft(:,end))) 1.1*max(weightft(:,end))])
% % title('Weighted summed spectrum')
% end

outfidw         =   weightsum(:,end);
outfidu         =   unweightsum(:,end);
outfid_wref     =   weightsum_wref(:,end);


