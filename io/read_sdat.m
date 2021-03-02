function [data,water,info,filename] = read_sdat()
% requires sdat in same folder
% get current directory
dircur = pwd;

% get the metabolite data
[filename,pathname,~]         = uigetfile('*_act.SDAT','Select .SDAT metabolite file');

cd(pathname);
% This is metabolite info
data	=   mrs_readSDAT(filename);
water	=   mrs_readSDAT([filename(1:end-8) 'ref.SDAT']); % contains water collected per dynamic
info    =   mrs_readSPAR([filename(1:end-5) '.SPAR']);

% convert SDAT files to spectral representation (and flip LR)
data    = mrs_fft(data);
water   = mrs_fft(water);

cd(dircur);

end

