function []=write_RAW_LCModel(path,fid,p_ID,header)

%####################################################
%               Save New File
%####################################################

%dlmwrite('old.RAW', metabo_spec2,'delimiter', '','precision','% 15.6e' );

W(:,1)=real(fid);
W(:,2)=imag(fid);


    
fileout = fopen(path,'w');
% fileout = fopen('/home/brian/.lcmodel/scripts_felix/LCModel_matlab_CSI/08d-11h-55m-10s-23008pid/met/RAW','w');

fprintf(fileout,' $SEQPAR\n');
fprintf(fileout,[' ECHOT = ',num2str(header.rdb_hdr.te/1e3),'\n']);%TE=68ms
fprintf(fileout,[' HZPPPM = ',num2str(header.rdb_hdr.ps_mps_freq/1e7,9),'\n']);% need accurate reading ftom header!
fprintf(fileout,' SEQ = ''PRESS''\n');
fprintf(fileout,' $END\n');
fprintf(fileout,' $NMID\n');
fprintf(fileout,' BRUKER = F,\n');
fprintf(fileout,[' ID = ''',p_ID,'''\n']);
fprintf(fileout,' FMTDAT = ''(2e15.6)'',\n');
fprintf(fileout,' VOLUME = 1.\n');%not relevant when using water scaling, needs adjusting if voxel volumes change
fprintf(fileout,' TRAMP = 1.\n');
fprintf(fileout,' BRUKER = F\n');
fprintf(fileout,' $END\n');
for index=1:size(fid,1)
    fprintf(fileout,'%15.6e%15.6e\n', W(index,1),W(index,2));
end

fprintf(fileout,'\n');

fclose('all');
end