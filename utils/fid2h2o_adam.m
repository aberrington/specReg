
function complexfid=fid2h2o_adam(directory, filename, vivofid)


fid(:,1)=real(vivofid);
fid(:,2)=imag(vivofid);


fileid=fopen([directory '/' filename '.H2O'],'w');
            fprintf(fileid,[' $NMID ID=''' filename ''', FMTDAT=''(8E13.5)''\n']);
            fprintf(fileid,[' TRAMP=1., VOLUME=1. $END\n']);
            fprintf(fileid,'%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n',[fid(:,1)'; fid(:,2)']); 
            fclose(fileid);

            
end