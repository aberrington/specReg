% File originally from KLChan Aug 2017 - modified by Adam Berrington
% Function to read in Philips data/list data and save them as a .spa file
% INPUTS:   
%           save_spa        output a .spa file, 1 or 0
%           fname           optional give a filename, or pick one from
%                           interface if none is given
% OUTPUTS:
%           metabfid        uncombined metabolite fid
%           waterfid        uncombined water fids
%           params          parameters for MRspa

function [metabfid, waterfid, params] = process_raw_philips(save_spa, fname, isHadamard)

    save_uncombined     =   1; % optional flag to save .mat data of uncombined coil data
    picksdat = 1;
    
if(nargin<3)
    isHadamard = 0;
end
    
if((nargin <2) | (fname == 0))
    dirs = uipickfiles('FilterSpec', '*.data', 'Prompt', 'Pick a .data file for reconstruction');
else
    dirs = {fname};
end

if(nargin<1)
    save_spa = 0;
end

if(save_spa || save_uncombined)
    % if saving a spa file - requires header from the corresponding SDAT
    % file to get relavant parameters

    a       = dir('*raw_act.SDAT');
    nfiles  = size(a,1);
    
    if(nfiles > 0 && picksdat ~= 1)
        disp(['Using 1st SDAT in folder for header information: ' a(1).name])
        dirs_sdat = {[a(1).folder '/' a(1).name]};
    else
        dirs_sdat = uipickfiles('Prompt', ['Please pick equivalent sdat file for information'], 'FilterSpec', '*.SDAT');
    end
end

for d = 1:length(dirs)

disp('Opening data...')

[a b]       = fileparts(dirs{d});
savename    = b;
this_file   = [b '.data'];


endian      = 'l';
type        = 'float';

disp('Reading scan parameters...')


fname_scan_params   = [this_file(1:(end-4)),'list']; % file name .list
scan_params         = textread(fname_scan_params, '%s'); % read whole .list file

% find the number of points in each FID
npoints_idx         = find(strcmp(scan_params,'t_range'));
npoints             = str2num(cell2mat(scan_params(npoints_idx(1)+3)))+1; % number of sampling FID points

% find the number of bytes for each FID
sizebyte_idx         = find(strcmp(scan_params,'size'));
sizebyte             = str2num(cell2mat(scan_params(sizebyte_idx(4)+(44)))); % number of sampling FID points

% find all instances of data and ignore the noise channels (NOI).
data_lines          = find(strcmp(scan_params,'STD')); % STD is standard data vector
tot_offset_idx      = scan_params(data_lines(end) + 20); % final STD + 20 more offsets - offset is the number of characters to get to the final offset in .data
data_lines          = data_lines(3:end); % ignore first 4 STDs - these are just descriptors
offset              = str2num(scan_params{data_lines(1) + 20}); % first data position
tot_offsets         = (str2num(tot_offset_idx{1}) - offset)/sizebyte + 1; % The number of offsets (position of the data) in .data (divided by 8192 bytes).

disp('Reading data...')

fp                  = fopen(this_file, 'rb', endian);
fseek(fp, offset, -1); % Move to the first offset of data

% read in the data (2048 * number of offsets)
data_raw            = fread(fp, npoints*2*tot_offsets, type); % npoints * 2 = number of complex points
                                       % size of each floating point. (unsigned
                                       % integer?
                                       
data                = data_raw(1:2:end,:)+1i*data_raw(2:2:end,:); % Points alternate between real and complex
data                = reshape(data,[npoints, tot_offsets]); % Reshape back to npoints by tot_offsets.

% number of lines of data (data_lines) = coils * (averages + water,spec)
D       = length(data_lines);

coil    = zeros(1,D-1);
kx      = coil;
ky      = coil;
avg     = coil;
sign    = coil;
on_off  = coil;
mix     = coil;

% extract info for each data line
nwater  = 0;

for dl = 1:D
    
    if dl == D
        this_dl = cellfun(@str2num,{scan_params{(data_lines(dl)+1):(data_lines(dl)+20)}});
    else
        this_dl = cellfun(@str2num,{scan_params{(data_lines(dl)+1):(data_lines(dl + 1)-1)}});
    end
        
    coil(dl)    = this_dl(6);   % record the coil number
    mix(dl)     = this_dl(1);   % record if water scan (for mixed sequences)
    avg(dl)     = this_dl(12);  % record the av. number
    on_off(dl)  = this_dl(7);   % record extra param 1 (=on/off for editing)
    sign(dl)    = this_dl(13);  % record sign of gradient direction
    
    if(mix(dl) == 1) % if it is a water scan
        nwater = max(nwater, avg(dl)); % store the number of water scans
    end

end

ncoils      = max(coil) + 1;
nav         = max(avg) + 1;
nwater      = nwater + 1;
non_off     = max(on_off);

on_count        = 0;
off_count       = 0;
all_dl          = [];
wat_ref_count   = 0;

if(non_off > 0)
    disp('there may be something different with the acquisiton - could be dynamic scans or edited spectra');
    metabfid = zeros(ncoils,nav, npoints, non_off);
    waterfid = zeros(ncoils,nav, npoints, non_off);
    
    for dl = 1:D
        is_water    = mix(dl);
        this_avg    = avg(dl)+1;
        this_coil   = coil(dl)+1;
        this_onoff  = on_off(dl)+1;

        if ~is_water % metabolite acquistion
            metabfid(this_coil, this_avg,:,this_onoff) = data(:,dl);
        else % water acquisiton
            waterfid(this_coil, this_avg,:,this_onoff) = data(:,dl);
        end
    end

else
    metabfid = zeros(ncoils, nav, npoints);
    waterfid = zeros(ncoils, nwater, npoints);

    for dl = 1:D
        is_water    = mix(dl);
        this_avg    = avg(dl)+1;
        this_coil   = coil(dl)+1;
        this_onoff  = on_off(dl)+1;

        if ~is_water % metabolite acquistion
            metabfid(this_coil, this_avg,:) = data(:,dl);
        else % water acquisiton
            waterfid(this_coil, this_avg,:) = data(:,dl);
        end
    end

end


out = io_loadspec_sdat(dirs_sdat{1},1);
params = [out.txfrq/1e6, out.spectralwidth, 4.65, 1, out.te, out.tr];
    
if(save_uncombined)
    save([savename '.mat'], 'metabfid', 'waterfid', 'params');
end

if(nav==1 && non_off>0) % if all dynamics no averages - treat as averages
    nav = non_off + 1;
    non_off = 0;
    metabfid = squeeze(permute(metabfid,[1,4,3,2]));
    waterfid = squeeze(permute(waterfid,[1,4,3,2]));

    nwaters_idx  =find(strcmp(scan_params,'number_of_extra_attribute_1_values'));
    nwater             = str2num(cell2mat(scan_params(nwaters_idx(2)+2)));
end

% take into account if there is dynamics - 
if(non_off > 0)
    M = metabfid;
    W = waterfid;
    clear metabfid waterfid
    
    outfid_ecc          = zeros(npoints, nav, non_off+1);
    outfid_wref_ecc     = zeros(npoints, nwater-1, non_off+1);
    
    for ndyn = 1:(non_off+1)
        
        metabfid = squeeze(M(:,:,:,ndyn));
        waterfid = squeeze(W(:,:,:,ndyn));

        outfid_weighted     = zeros(npoints, nav);
        outfid_unweighted   = zeros(npoints, nav);

        disp('Weighting and ECC average number... ');

        for n = 1:nav
            disp(num2str(n));
            [outfid_weighted(:,n), outfid_unweighted(:,n), weights, outfid_wref(:,1)]   =   sum_phased_array(squeeze(metabfid(:,n,:))', squeeze(waterfid(:,1,:))', 1, 0);
            outfid_ecc(:,n,ndyn) = spec_ecc_cor(outfid_weighted(:,n), outfid_wref(:,1),0,0);
        end

        disp('Weighting and ECC average water scan... ');
        for m = 2:nwater
            % correct remaining water scans
            disp(num2str(m-1));
            [outfid_wref(:,m), ~, ~, ~]         =   sum_phased_array(squeeze(waterfid(:,m,:))', squeeze(waterfid(:,1,:))', 1, 0);
            outfid_wref_ecc(:,m-1,ndyn)              =   spec_ecc_cor(outfid_wref(:,m),outfid_wref(:,1),0,0); % have to now throw away 1st water scan
        end

    end
else
    outfid_weighted     = zeros(npoints, nav);
    outfid_unweighted   = zeros(npoints, nav);
    outfid_ecc          = zeros(npoints, nav);
    outfid_wref_ecc     = zeros(npoints, nwater-1);
    disp('Weighting and ECC average number... ');

    for n = 1:nav
        disp(num2str(n));
        [outfid_weighted(:,n), outfid_unweighted(:,n), weights, outfid_wref(:,1)]   =   sum_phased_array(squeeze(metabfid(:,n,:))', squeeze(waterfid(:,1,:))', 1, 0);
        outfid_ecc(:,n) = spec_ecc_cor(outfid_weighted(:,n), outfid_wref(:,1),0,0);
    end

    disp('Weighting and ECC average water scan... ');
    
    if(nwater==1)
        disp('Only 1 water reference found')
        [outfid_wref(:,1), ~, ~, ~]         =   sum_phased_array(squeeze(waterfid(:,1,:))', squeeze(waterfid(:,1,:))', 1, 0);
        outfid_wref_ecc(:,1)                =   spec_ecc_cor(outfid_wref(:,1),outfid_wref(:,1),0,0); % have to now throw away 1st water scan
    end
    
    for m = 2:nwater
        % correct remaining water scans
        disp(num2str(m-1));
        [outfid_wref(:,m), ~, ~, ~]         =   sum_phased_array(squeeze(waterfid(:,m,:))', squeeze(waterfid(:,1,:))', 1, 0);
        outfid_wref_ecc(:,m-1)              =   spec_ecc_cor(outfid_wref(:,m),outfid_wref(:,1),0,0); % have to now throw away 1st water scan
    end
    
    % if Hadamard then need to decode and do ecc
    
    if(isHadamard)
        outfid_wref_ecc_v1 = [];
        outfid_wref_ecc_v2 = [];
        
        metabv1 = outfid_weighted(:,1:2:end) + outfid_weighted(:,2:2:end);
        metabv2 = outfid_weighted(:,1:2:end) - outfid_weighted(:,2:2:end);
        
        waterv1 = outfid_wref(:,1:2:end) + outfid_wref(:,2:2:end);
        waterv2 = outfid_wref(:,1:2:end) - outfid_wref(:,2:2:end);
        
        for n = 1:(nav/2)
            outfid_ecc_v1(:,n) = spec_ecc_cor(metabv1(:,n), waterv1(:,1),0,0);
            outfid_ecc_v2(:,n) = spec_ecc_cor(metabv2(:,n), waterv2(:,1),0,0);
        end
        
        for m = 2:(nwater/2)
            outfid_wref_ecc_v1(:,m-1) = spec_ecc_cor(waterv1(:,m),waterv1(:,1),0,0); % have to now throw away 1st water scan
            outfid_wref_ecc_v2(:,m-1) = spec_ecc_cor(waterv2(:,m),waterv2(:,1),0,0); % have to now throw away 1st water scan
        end
    end
    
end

metabfid    = outfid_ecc;
waterfid    = outfid_wref_ecc;
params      = [out.txfrq/1e6, out.spectralwidth, 4.65, 1, out.te, out.tr];

%% save as spa file
if(save_spa)
clear data

%dirs_sdat = uipickfiles('Prompt', ['Please pick equivalent sdat file for ' b ' information'], 'FilterSpec', '*.SDAT');
out = io_loadspec_sdat(dirs_sdat{1},1);
disp(['Spectral width found to be:' out.spectralwidth]);
data.metab              = outfid_ecc;
data.arrayedmetab       = outfid_ecc;
data.ws                 = outfid_wref_ecc; % sum(outfid_wref_ecc,2);
data.ntmetab            = nav;
data.ntws               = size(outfid_wref_ecc,2);
data.params             = [out.txfrq/1e6, out.spectralwidth, 4.65, 1]; % taken from .sdat but need to write in future?
data.waterTEs.CSFfrac   = 0; 

save([savename '_ecc.spa'], 'data');

data.metab              = outfid_weighted;
data.arrayedmetab       = outfid_weighted;

save([savename '.spa'], 'data');


if(isHadamard)
    data.metab              = outfid_ecc_v1;
    data.arrayedmetab       = outfid_ecc_v1;
    data.ws                 = outfid_wref_ecc_v1;
    data.ntmetab            = nav/2;
    save([savename '_ecc_v1.spa'], 'data');
    
    data.metab              = outfid_ecc_v2;
    data.arrayedmetab       = outfid_ecc_v2;
    data.ws                 = outfid_wref_ecc_v2;
    data.ntmetab            = nav/2;
    save([savename '_ecc_v2.spa'], 'data');
end

end


end

