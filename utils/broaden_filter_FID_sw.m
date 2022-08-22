function [fid_b]=broaden_filter_FID_sw(fid, LB_factor, sw, filter_factor)
%Does line broadening and also filtering as gaussian
% Input: fid                 FID 
%        LB_factor           factor in Hz
%        filter_factor       Gaussian Filter 
%
% Output:fid_b               broadened FID


%sw  =   6000; % Assume SW = 6000Hz

if(isnan(LB_factor))
    LB_factor = 0;
end
dw  =   1/sw;
t   =   [0:dw:dw*(length(fid)-1)]';

fid_b = zeros(size(fid));

if(nargin < 4)
    for i = 1:size(fid,2)
        fid_b(:,i) = fid(:,i).*exp(-t*pi*LB_factor);
    end
else
    for i = 1:size(fid,2)
        fid_b(:,i) = fid(:,i).*exp(-t*pi*LB_factor - t.^2/(filter_factor^2));
    end
end
