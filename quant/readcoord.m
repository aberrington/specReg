function [lcmoutput, nummetab] = readcoord(filename)

%nummetab=24;
fid = -1;
fid=fopen(filename,'r');

if(fid < 0)
    for j = 1:50
        lcmoutput.metabconc(j).absconc = NaN;
        lcmoutput.metabconc(j).relconc = NaN;
        lcmoutput.metabconc(j).SD = NaN;
        lcmoutput.metabconc(j).name = NaN;
    end
    lcmoutput.linewidth = NaN;
    lcmoutput.SN = NaN;
else
    
    nummetab=textscan(fid,'%f',1,'EmptyValue',0,'headerlines',2);
    
    nummetab = cell2mat(nummetab);
    nummetab = nummetab -1;
    
    disp(['Number of metabolite fits found ' num2str(nummetab)])
    
    data=textscan(fid,'%f %s %f %s',nummetab,'EmptyValue',0,'headerlines',2);
    
    AB=data{1}; %absolute conc.
    SD=data{2}; %CRLB
    RL=data{3}; %relative conc.
    NM=data{4}; %metabolite names
    
    if(size(SD,1) == 1 || size(SD,1) == 0)
        
        for j=1:nummetab
            lcmoutput.metabconc(j).absconc = NaN;
            lcmoutput.metabconc(j).relconc = NaN;
            lcmoutput.metabconc(j).SD = NaN;
            lcmoutput.metabconc(j).name = NaN;
        end
        
        lcmoutput.linewidth = NaN;
        lcmoutput.SN = NaN;
        
    else
        SD_temp=cell(nummetab,1);
        
        for i=1:nummetab
            SD_temp{i}=sscanf(SD{i},['%f', char(176)]);
        end
        
        SD=cell2mat(SD_temp);
        
        for j=1:nummetab
            
            lcmoutput.metabconc(j).absconc = AB(j);
            lcmoutput.metabconc(j).relconc = RL(j);
            lcmoutput.metabconc(j).SD = SD(j);
            lcmoutput.metabconc(j).name = NM(j);
        end
        
        data=textscan(fid,'%s %s %f %s %s %s %f', 1,'EmptyValue',0,'headerlines',2);
        
        lcmoutput.linewidth = data{3};
        lcmoutput.SN = data{7};
        
    end
    fclose(fid);
end


  
