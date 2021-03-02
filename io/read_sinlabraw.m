function [data,water,info,filename] = read_sinlabraw()
% requires sdat in same folder
% get current directory
dircur = pwd;

[filename1,~,~]         = uigetfile('.raw','Select .raw file');
[filename,pathname,~]  = uigetfile('.SPAR','Select .SPAR file');

cd(pathname);

[data, water]       =   mrs_readRAW(filename1); % reads in the sin/lab/raw data and gets the metab/water specs
info                =   mrs_readSPAR(filename);

cd(dircur);

end

