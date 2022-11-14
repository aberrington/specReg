function h = read_header_ox(fname)
%READ_HEADER_OX  Read Pfile or ScanArchive header via Orchestra
%         h = read_header_value(fname)
%     fname   filename
%         h   Header structure
%
% 1/2021  Rolf Schulte
if (nargout<1), help(mfilename); return; end


%% check for Ox
if ~exist('GERecon','file')
    warning('GERecon not found');
end


%% check if file exists
if ~exist(fname,'file')
    warning('fname (=%s) not found',fname);
end


%% load header
if ~isempty(regexpi(fname,'\.7$'))
    pfile = GERecon('Pfile.Load', fname);
    h = GERecon('Pfile.Header',pfile);
    h = pfile2header(h);
else
    if ~isempty(regexpi(fname,'\.h5$'))
        archive = GERecon('Archive.Load', fname);
        GERecon('Archive.Close', archive);
        h = archive2header(archive);
    else
        error('fname (=%s) not ending with .7 or .h5',fname);
    end 
end


end      % read_header_ox.m
