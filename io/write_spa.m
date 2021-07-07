function data = write_spa(pathname, fid, waterfid, txfrq, sw, delta0) 

data.metab      = fid;%(out.spectralaverages)*out.fids'; % scale it to match .spa data since sdat is average NOT sum
data.ws         = waterfid;
data.ntmetab    = 1; % do everything 'averaged' so this =1;
data.ntws       = 1; % do everything 'averaged' so this =1;
data.params(1)  = txfrq / 1e6;
data.params(2)  = sw;
data.params(3)  = delta0; % usually 4.65
data.params(4)  = 1;
data.sdatinfo   = [];
data.waterTEs.CSFfrac = 0;

save(pathname, 'data');

end
