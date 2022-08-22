% Script to run LCModel for GABA data

% want to find the .spa files in the SpecReg folder and then run them
% through LCModel

% Make an LCModel folder in SpecReg
filepath = 'SpecReg/LCModel';
mkdir(filepath);

% Get all the processed files
D       =   dir('SpecReg/Data/*.spa');
C       =   strsplit(D(1).name,'_');

b = dir('SpecReg/Data/Block/'); % Search for block (fMRS analysis)
numBlockAnalyses = length(b)-2;

if(numBlockAnalyses>0) % 3rd position since .. and . are directories...
    for j = 1:numBlockAnalyses
        BlockName{j} = b(j+2).name;
    end
else
    numBlockAnalyses =0;
end

fName{1} = ['SpecReg/Data/' C{1}];
fPath{1} = 'SpecReg/LCModel';
    
for k = 1:(numBlockAnalyses)
    fName{k+1} = ['SpecReg/Data/Block/' BlockName{k} '/' C{1}];
    fPath{k+1} = ['SpecReg/LCModel/Block/' BlockName{k}];
    mkdir(fPath{k+1});
end

specType = {'diff', 'off'};

for n = 1:(numBlockAnalyses + 1)
    
    for i = 1:2

        fileName = [fName{n} '_' specType{i}];

        load([fileName '.spa'], '-mat');

        nSpec = size(data.metab,2);
        
        for s = 1:nSpec
           
            if(nSpec == 1)
                savename = [C{1} '_' specType{i}];
            else
                savename = [C{1} '_' specType{i} '_' num2str(s)];
            end

        specData    = data.metab(:,s); % spec data
        waterData   = data.ws; % water data
        nbpoints    =   size(specData, 1); % number of samples
        % these should = 1
        ntmetab     =   data.ntmetab; % only 1 metab - seems to be averaged to 1: specData.spectralaverages; % number of averages
        ntwater     =   data.ntws; % need to check this - not accessible from SPAR on Examcard (Spectral Corr NSA)

        sw          =   data.params(2);
        sfrq        =   abs(data.params(1)); % in Hz / ppm for LCModel
        
        if(data.params(1) > 290) % This is 7T
           fieldS = 7;
        else
           fieldS = 3; 
        end
        
        if(fieldS==3)
            TE          =   68; % 7T is 72, 3T is 68
            specData    =   specData'; % currently need transpose for 3T GE dataset
            waterData   =   waterData';
        else
            TE          =   72;
        end
        
        fid2raw_adam(fPath{n}, savename, specData') % IDX may change
        fid2h2o_adam(fPath{n}, savename, waterData');

        % Change to where basis is located on server
        if(i==1)
            if(fieldS==3)
                tes_basis(1).te     = '/home/bbzapb/bases/3t_GE_MEGAPRESS_june2011_diff';
            else
               tes_basis(1).te     = '/home/bbzapb/bases/7T_MEGAsLASER_72ms_DIFF_2021';
            end
        else
            if(fieldS==3)
                tes_basis(1).te     = '/home/bbzapb/bases/3T_PRESS68_no_LMM';
            else
                tes_basis(1).te     = '/home/bbzapb/bases/7T_MEGAsLASER_72ms_OFF_2021';
            end
        end
        
        % Concentration corrections - these are not done on GE data... leave
        % them here
%         T2metab     =   9999999; %values from Lactate at 7T
%         T2water     =   88; %values from Michaeli
%         tissuewatercontent  = 0.819;  % tissue water content may vary depending on the brain region. (eg 0.71 for pure white matter, 0.81 for pure grey matter)
%         CSFwatercontent     = 1;      % assuming CSF is pure water
%         CSFfrac             = 0;      % CSF fraction
%         attmet              = ntmetab*exp(-TE/T2metab);
%         concpurewater       = 55555;  % concentration of pure water in mM
%         atth2o              = ntwater*exp(-TE/T2water);
%         atth2o              = atth2o(:,1);
%     %     %factor of 1.15 from Michaeli MRM 2003 estimated at 3T
%     %     %atth2o = data.ntws;%comment
%         wconc   = concpurewater * (CSFfrac*CSFwatercontent + (1-CSFfrac)*tissuewatercontent);
%         wconc   = wconc./(1-CSFfrac);

        % Make LCModel Control file
        fileid      =   fopen([fPath{n} '/' savename '.CONTROL'],'w');

        fprintf(fileid,' $LCMODL\n');
        fprintf(fileid,[' TITLE=''' savename '''\n']);
        fprintf(fileid,' KEY=210387309\n');  % Key for Bernoulli
        fprintf(fileid,' PGNORM=''US''\n');
        fprintf(fileid,[' FILPS=''' savename '.PS''\n']);
        fprintf(fileid,[' FILCOO=''' savename '.COORD''\n']);
        fprintf(fileid,[' ECHOT=' num2str(TE) '\n']);
        
        if(i==1)
            print_control(fileid, fieldS, 'MEGA_DIFF');
        else
            print_control(fileid, fieldS, 'MEGA_OFF');
        end
        
        
        %fprintf(fileid,' IAVERG=3\n');
        %fprintf(fileid,' VITRO=T\n');
%         fprintf(fileid,' PPMSHF=0.0\n');
%         fprintf(fileid,' PPMEND=0.2\n');
%         fprintf(fileid,' PPMST=4.0\n');
%         if(i==1)
%         fprintf(fileid,' PPMGAP(1,1)=1.95\n');
%         fprintf(fileid,' PPMGAP(2,1)=1.2\n');
%         end
%         
%         if(i==1)
%             %fprintf(fileid,' DKNTMN=999\n');
%         end
%             
%         fprintf(fileid,' RFWHM=3\n');
%         fprintf(fileid,' FWHMBA=0.0050\n');
%         fprintf(fileid,' SHIFMX=0.3,0.3\n');
%         fprintf(fileid,' SHIFMN=-0.2,-0.1\n');
%         fprintf(fileid,' SDDEGP=5.00\n');
%         fprintf(fileid,' SDDEGZ=999.00\n');
%         fprintf(fileid,' NEACH=99\n');
%         %fprintf(fileid,' CHOMIT(6)=''bHB''\n');
%         %fprintf(fileid,' CHOMIT(5)=''Cho''\n');
%         %fprintf(fileid,' CHOMIT(4)=''Thr''\n');
%         %fprintf(fileid,' CHOMIT(3)=''Gly''\n');
%         %fprintf(fileid,' CHOMIT(2)=''ATP''\n');
%         %fprintf(fileid,' CHOMIT(1)=''Ace''\n');
%         %fprintf(fileid,' NOMIT=6\n');
%         fprintf(fileid,' LCOORD=9\n');
%         fprintf(fileid,' LPRINT=6\n');
%         %fprintf(fileid,' NRATIO=0\n');
%         fprintf(fileid,' DEGZER=0.0\n');
%         fprintf(fileid,' NCOMBI=3\n');
%         fprintf(fileid,' CHCOMB(1)=''Glu+Gln''\n');
%         fprintf(fileid,' CHCOMB(2)=''Glu+Gln+GSH''\n');
%         fprintf(fileid,' CHCOMB(3)=''NAA+NAAG''\n');
%         %fprintf(fileid,' CHCOMB(4)=''Glu+Gln''\n');
%         %fprintf(fileid,' CHCOMB(5)=''Glc+Tau''\n');
%         fprintf(fileid,' CONREL=1.00\n');
%         fprintf(fileid,' NAMREL=''NAA+NAAG''\n');
%         %fprintf(fileid,' CHUSE1(1)=''Cr''\n');
%         %fprintf(fileid,' CHUSE1(2)=''PCr''\n');
%         %fprintf(fileid,' CHUSE1(3)=''NAA''\n');
%         %fprintf(fileid,' CHUSE1(4)=''Glu''\n');
%         %fprintf(fileid,' CHUSE1(5)=''Ins''\n');
%         %fprintf(fileid,' NUSE1=5\n');
%         fprintf(fileid,' DEGPPM=0.\n');
%         fprintf(fileid,' DEGZER=0.\n');
%         %fprintf(fileid,' NSIMUL=0\n');
%         %fprintf(fileid,' NCALIB=0\n');
%         %fprintf(fileid,' WSPPM=3.0241\n');
%         %fprintf(fileid,' WSPPM=2.01\n');
%         %fprintf(fileid,' WSMET=''NAA''\n');
%         %fprintf(fileid,' N1HMET=3\n');
%         
% 
%         % not in STANDARD.CONTROL
%         %fprintf(fileid,' NSIMUL=0\n');
%         fprintf(fileid,' DOWS=T\n');
%         fprintf(fileid,' DOECC=F\n');
% 
%         % using default values for corrections
%         %fprintf(fileid,[' ATTH2O=' num2str(atth2o*1) '\n']);
%         %fprintf(fileid,[' ATTMET=' num2str(2) '\n']); % 2 when difference spectrum
%         fprintf(fileid,[' ATTH2O=' num2str(0.7) '\n']); % to match current output
%         fprintf(fileid,[' WCONC=' num2str(35880) '\n']);
        %fprintf(fileid,[' WCONC=' num2str(wconc) '\n']);
        fprintf(fileid,[' FILPRI=''' savename '.PRINT''\n']);
        fprintf(fileid,[' FILRAW=''' savename '.RAW''\n']);
        fprintf(fileid,[' FILH2O=''' savename '.H2O''\n']);
        fprintf(fileid,[' FILBAS=''' tes_basis(1).te '.basis''\n']);
        fprintf(fileid,[' DELTAT=' num2str(1/sw) '\n']);
        fprintf(fileid,[' NUNFIL=' num2str(nbpoints) '\n']);%%nbpoints should be defined
        fprintf(fileid,[' HZPPPM=' num2str(sfrq) '\n']);
        fprintf(fileid,' $END\n');
        fclose(fileid);

        % run LCModel
        end

    end

end

