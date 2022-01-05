

    % read in the COORD file
    lcmoutput    = readcoord_adam([filepath '/' savename '.COORD']);

    metabidx = [];
    for mm = 1:length(lcmoutput.metabconc)
        if(strcmp(lcmoutput.metabconc(mm).name,'GABA'))
            metabidx = mm;
        end
    end
    
    if(~isempty(metabidx))
        % choose ABSOLUTE conc
        conc = lcmoutput(n).metabconc(metabidx).absconc;
    else
        conc = NaN;
    end
        