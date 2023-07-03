% function to read in data list
function [data water info filename] = read_datalist()

if nargin==0
    [filename,filepath]=uigetfile('*.data', 'Select .data file');
end

[data_fid, water_fid, params] = process_raw_philips(0, filename);

% convert to spectra
data    = mrs_fft(conj(data_fid));
water   = mrs_fft(conj(water_fid));

info.BW = params(2);
info.transmit_frequency = params(1)*10^6;
info.TE = params(5);
info.TR = params(6);

end