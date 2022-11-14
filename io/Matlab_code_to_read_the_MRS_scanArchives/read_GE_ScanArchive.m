% Function to read ScanArchive data and output same information as P-file

function [data, water, info, filename] = read_GE_ScanArchive()

zf_factor=1;

%% read data
if nargin==0
    [filename,filepath]=uigetfile('*.h5', 'Select GE ScanArchive');
end

[d,h,archive]=read_archive([filepath, filename]);

disp('-----------')
disp([filename,':'])
get_MRS_info(h);

all_spectra = d;


end
