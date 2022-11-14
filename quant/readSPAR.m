function [TE, TR, B0] = readSPAR(~)

%nummetab=24;
fid = -1;
dirs = uipickfiles('Prompt', 'Choose corresponding MEGA SPAR file', 'FilterSpec', '*.SPAR');

s=fileread(dirs{1});
TE = regexp(s,'echo_time : (\w*)','tokens');
TE = str2num(TE{1}{1});
TR = regexp(s,'repetition_time : (\w*)','tokens');
TR = str2num(TR{1}{1});
f0 = regexp(s,'synthesizer_frequency : (\w*)','tokens');
f0 = str2num(f0{1}{1});
B0 = round(f0/42.577*1e-6);





