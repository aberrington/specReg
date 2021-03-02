% Script to run LCModel for GABA data

% want to find the .spa files in the SpecReg folder and then run them
% through LCModel

% Make an LCModel folder in SpecReg
mkdir('SpecReg/LCModel');

% Get all the processed files
D   =   dir('SpecReg/Data/*.spa');

load(dirs{i}, '-mat');

[a, b] = fileparts(dirs{i});
savename = b;
%fid_b=broaden_filter_FID_sw(data.metab, 4, data.params(2));
specData    = data.metab; % spec data
waterData   = data.ws; % water data

%% Define parameters
nbpoints    =   size(specData, 1); % number of samples
ntmetab     =   data.ntmetab; % only 1 metab - seems to be averaged to 1: specData.spectralaverages; % number of averages
ntwater     =   data.ntws; % need to check this - not accessible from SPAR on Examcard (Spectral Corr NSA)
                   % seems like you don't need this with
                   % Philips processed files. Maybe the average
                   % is just the arithmetic mean not simply
                   % summing but should be 4 in this study
sw          =   data.params(2);
sfrq        =   data.params(1); % in Hz / ppm for LCModel
TE          =   144;

T2metab     =   9999999; %values from Lactate at 7T
T2water     =   88; %values from Michaeli

tissuewatercontent  = 0.819;  % tissue water content may vary depending on the brain region. (eg 0.71 for pure white matter, 0.81 for pure grey matter)
CSFwatercontent     = 1;      % assuming CSF is pure water
CSFfrac             = 0;      % CSF fraction

attmet              = ntmetab*exp(-TE/T2metab);

concpurewater       = 55555;  % concentration of pure water in mM
atth2o              = ntwater*exp(-TE/T2water);
atth2o              = atth2o(:,1);

%factor of 1.15 from Michaeli MRM 2003 estimated at 3T
%atth2o = data.ntws;%comment

wconc   = concpurewater * (CSFfrac*CSFwatercontent + (1-CSFfrac)*tissuewatercontent);
wconc   = wconc./(1-CSFfrac);

%plot(real(fftshift(fft(conj(specData.fids)))))
fid2raw_adam(filepath, savename, specData) % IDX may change
fid2h2o_adam(filepath, savename, waterData);

%tes_basis(1).te     = '/home/bbzapb/bases/STEAM_TE16_TM16_7T'; % directory on radpub
%tes_basis(1).te     = '/home/bbzapb/bases/7T/slaser_38ms/real/RAW/7T_sLASER_38ms_REAL_2019';
tes_basis(1).te     = '/home/bbzapb/bases/sLASER_MM';
% newdir=[filename '_lcm' ];
% mkdir(directory,newdir);
% lcmodeloutputfilename=['lcm_' num2str(i)PGNORM = 'US'

fileid      =   fopen([filepath '/' savename '.CONTROL'],'w');

fprintf(fileid,' $LCMODL\n');
fprintf(fileid,[' TITLE=''' savename '''\n']);
%fprintf(fileid,' OWNER=''Johns Hopkins University, Neuroradiology''\n');
%fprintf(fileid,' KEY=217638264\n');  % Key for Bernoulli
fprintf(fileid,' PGNORM=''US''\n');
fprintf(fileid,[' FILPS=''' savename '.PS''\n']);
fprintf(fileid,[' FILCOO=''' savename '.COORD''\n']);

fprintf(fileid,' PPMSHF=0.0\n');
%fprintf(fileid,' PPMEND=0.5\n');
%fprintf(fileid,' PPMST=4.2\n');
% for editing

fprintf(fileid,' PPMEND=0.5\n');
fprintf(fileid,' PPMST=4.2\n');

fprintf(fileid,' DKNTMN=99\n');
%fprintf(fileid,' RFWHM=2.5\n');
fprintf(fileid,' FWHMBA=0.0050\n');
fprintf(fileid,' SHIFMX=0.3,0.3\n');
fprintf(fileid,' SHIFMN=-0.2,-0.1\n');
fprintf(fileid,' SDDEGP=5.00\n');
fprintf(fileid,' SDDEGZ=5.00\n');
fprintf(fileid,' NEACH=99\n');

fprintf(fileid,' CHOMIT(6)=''bHB''\n'); 
fprintf(fileid,' CHOMIT(5)=''Cho''\n');  
fprintf(fileid,' CHOMIT(4)=''Thr''\n'); 
fprintf(fileid,' CHOMIT(3)=''Gly''\n');
fprintf(fileid,' CHOMIT(2)=''ATP''\n');
fprintf(fileid,' CHOMIT(1)=''Ace''\n');
fprintf(fileid,' NOMIT=6\n');

fprintf(fileid,' LCOORD=9\n');
fprintf(fileid,' LPRINT=6\n');
fprintf(fileid,' NRATIO=0\n'); 

fprintf(fileid,' DEGZER=0.0\n');
fprintf(fileid,' NCOMBI=5\n');
fprintf(fileid,' CHCOMB(1)=''PCho+GPC''\n');
fprintf(fileid,' CHCOMB(2)=''Cr_1+Cr_2+PCr_1+PCr_2''\n');
fprintf(fileid,' CHCOMB(3)=''NAA_s+NAA_m+NAAG''\n');
fprintf(fileid,' CHCOMB(4)=''Glu+Gln''\n');
fprintf(fileid,' CHCOMB(5)=''Glc+Tau''\n');

fprintf(fileid,' CONREL=1.00\n');
fprintf(fileid,' NAMREL=''Cr_1+Cr_2+PCr_1+PCr_2''\n');
fprintf(fileid,' CHUSE1(1)=''Cr''\n');
fprintf(fileid,' CHUSE1(2)=''PCr''\n');
fprintf(fileid,' CHUSE1(3)=''NAA''\n');
fprintf(fileid,' CHUSE1(4)=''Glu''\n');
fprintf(fileid,' CHUSE1(5)=''Ins''\n');
fprintf(fileid,' NUSE1=5\n');
fprintf(fileid,' DEGPPM=0.\n');
fprintf(fileid,' DEGZER=0.\n');
fprintf(fileid,' NSIMUL=0\n');

%fprintf(fileid,' NCALIB=0\n');
fprintf(fileid,' WSPPM=3.0241\n');
fprintf(fileid,' WSMET=''Cr_2''\n');
fprintf(fileid,' N1HMET=3\n');

fprintf(fileid,[' FILPRI=''' savename '.PRINT''\n']);

% not in STANDARD.CONTROL
%fprintf(fileid,' NSIMUL=0\n');

fprintf(fileid,' DOWS=T\n');
fprintf(fileid,' DOECC=F\n');
fprintf(fileid,[' ATTH2O=' num2str(atth2o*1) '\n']);
fprintf(fileid,[' ATTMET=' num2str(attmet) '\n']);
fprintf(fileid,[' WCONC=' num2str(wconc) '\n']);
%%%%%%
fprintf(fileid,[' FILRAW=''' savename '.RAW''\n']);
fprintf(fileid,[' FILH2O=''' savename '.H2O''\n']);
fprintf(fileid,[' FILBAS=''' tes_basis(1).te '.basis''\n']);
fprintf(fileid,[' DELTAT=' num2str(1/sw) '\n']);
fprintf(fileid,[' NUNFIL=' num2str(nbpoints) '\n']);%%nbpoints should be defined
fprintf(fileid,[' HZPPPM=' num2str(sfrq) '\n']);
fprintf(fileid,' $END\n');
fclose(fileid);
