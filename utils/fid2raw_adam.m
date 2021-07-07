
function complexfid=fid2raw_adam(directory, filename, vivofid)

fid(:,1)=real(vivofid);
fid(:,2)=imag(vivofid);

fileid=fopen([directory '/' filename '.RAW'],'w');
            fprintf(fileid,[' $NMID ID=''' filename ''', FMTDAT=''(8E13.5)''\n']);
            fprintf(fileid,[' TRAMP=1., VOLUME=1. $END\n']);
            %fprintf(fileid,'    %#.6E    %#.6E\n',[fid(:,1)'; fid(:,2)']); 
            fprintf(fileid,'%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n',[fid(:,1)'; fid(:,2)']); 
            fclose(fileid);

            
end