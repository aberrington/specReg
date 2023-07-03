% script will run tissue segmentation using equations in 
%% https://onlinelibrary.wiley.com/doi/10.1002/nbm.4257
function [] = specreg_quant_STEAM(tissue_frac)

if(nargin<1)
    disp('tissue fractions not provided - will look for .json file')
    flag_frac = 0;
else
    flag_frac = 1;
    disp(['GM is ' num2str(tissue_frac(1)*100) '%'])
    frac.GM = tissue_frac(1);
    disp(['WM is ' num2str(tissue_frac(2)*100) '%'])
    frac.WM = tissue_frac(2);
    disp(['CSF is ' num2str(tissue_frac(3)*100) '%'])
    frac.CSF = tissue_frac(3);
end

filepath = 'SpecReg/STEAM/LCModel';

% fit diff and also off
%% Get sequence parameters
infofileID = 'SpecReg/STEAM/info.json';
val = read_json_info(infofileID);

TE = val.TE;
TR = val.TR;
B0 = val.B0;

disp(['Found TE/TR (ms): ' num2str(TE) '/' num2str(TR)]);

mkdir('SpecReg/STEAM/Quant');
    fileName = ['raw'];

    [lcmoutput, nmetabs] = readcoord([filepath '/' fileName '.COORD']);


    % print csv file
    %% Read water T2 values - expected in the working folder
    fid = fopen('tissue_values.json'); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    density = jsondecode(str);
    %% Water relaxation - use default 7 T values
    disp(['Using default ' num2str(B0) 'T values'])
    fid = fopen(['RH2O_' num2str(B0) 'T.json']);
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    relaxH2O = jsondecode(str);

    disp(['Using H20 T1 times (GM, WM, CSF) in ms: ' ...
        num2str(relaxH2O.T_GM(1)) ', ' num2str(relaxH2O.T_WM(1)) ', ' num2str(relaxH2O.T_CSF(1))]);
    disp(['Using H20 T2 times (GM, WM, CSF) in ms: ' ...
        num2str(relaxH2O.T_GM(2)) ', ' num2str(relaxH2O.T_WM(2)) ', ' num2str(relaxH2O.T_CSF(2))]);
    %% Metabolite relaxation
    fid = fopen(['RM_' num2str(B0) 'T.json']);
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    relaxM = jsondecode(str);

    %% Tissue fraction
    if(flag_frac == 0)
        fid = fopen('tissue_frac_segmentation.json');
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        frac = jsondecode(str);
        disp(['Using tissue fractions (GM, WM, CSF): ' ...
            num2str(frac.GM(1)) ', ' num2str(frac.WM(1)) ', ' num2str(frac.CSF(1))]);
    end


    %% Calculate relaxation coefficients
    RH2O.WM = exp(-TE/relaxH2O.T_WM(2)) * (1-exp(-TR/relaxH2O.T_WM(1)));
    RH2O.GM = exp(-TE/relaxH2O.T_GM(2)) * (1-exp(-TR/relaxH2O.T_GM(1)));
    RH2O.CSF = exp(-TE/relaxH2O.T_CSF(2)) * (1-exp(-TR/relaxH2O.T_CSF(1)));

    CRLB        = [lcmoutput.metabconc.SD]';
    Metabolite  = convertCharsToStrings([lcmoutput.metabconc.name])';
    LCM_Relative     = [lcmoutput.metabconc.relconc]';
    LCM_Conc    = [lcmoutput.metabconc.absconc]';


    % Haven't included metabolite relaxation
    for n = 1:nmetabs

        if(strcmp(Metabolite(n),"NAA+NAAG"))
            m_string = "tNAA";
        elseif(strcmp(Metabolite(n), "Glu+Gln"))
            m_string = "Glx";
        elseif(strcmp(Metabolite(n), "Cr+PCr"))
            m_string = "tCr";
        elseif(strcmp(Metabolite(n), "PCho+GPC"))
            m_string = "tCho";
        else
            m_string = Metabolite(n);
        end

        RMETAB = calc_relax_metab(m_string, relaxM, TE, TR);

        Scaled_Conc(n) = LCM_Conc(n) * (frac.GM(1)*density.d_GM*RH2O.GM ... 
                + frac.WM(1)*density.d_WM*RH2O.WM ...
                + frac.CSF(1)*density.d_CSF*RH2O.CSF)/((1-frac.CSF)*RMETAB);
        
        Metabolite(n) = m_string;
        
        if(RMETAB == -1)
            EstConc{n} = ' ';
        else
            EstConc{n} = Scaled_Conc(n);
        end
    end


    
    varNames = {'Metabolite','CRLB','Relative Concentration', 'LCModel Concetration', 'Estimated Tissue Corrected Concentration (mM)'};
    T = table(Metabolite, CRLB, LCM_Relative, LCM_Conc, EstConc', 'VariableNames', varNames);
    writetable(T, ['SpecReg/STEAM/Quant/concs.csv'])

    % read in file raw_diff and raw_off
    % perform corrections
end

% write to a csv file
function RMETAB = calc_relax_metab(metabname, relaxM, TE, TR)

    if(~isfield(relaxM, metabname))
        disp(['No relaxation values for ' metabname])
        RMETAB = -1;
    else
        RMETAB = exp(-TE/relaxM.(metabname)(2)) * (1-exp(-TR/relaxM.(metabname)(1)));
    end
end
