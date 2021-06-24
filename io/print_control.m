function print_control(fileid, fieldS, seq)

control_file = [num2str(fieldS) 'T_' seq '.txt'];
fid = fopen(control_file, 'r');
tline = fgetl(fid);
while tline ~= -1
    fprintf(fileid, [tline '\n']);
    tline = fgetl(fid);
end

fclose(fid);
end
