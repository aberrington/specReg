function fid = eddyCor2(fidw,fidm)
%
% function fid=eddyCor2(waterfid,metabfid)
% DC, ECC and zero-order phase corrections
% 
% Dinesh Deelchand 21 March 2008
%
% Updated for array experiment
% Added DC correction after ECC - 03 Nov 2011
% Updated for array experiment
% Updated now removes phase of first point (based on Ivan Tkac macro) - 21 March 2012

if (size(fidw,2)==1)
    [~,nt] = size(fidm);
    fidw = repmat(fidw,[1 nt]);
end

% DC correction
fidw = fidw - mean(fidw(length(fidw)-200:length(fidw)));
fidm = fidm - mean(fidm(length(fidm)-200:length(fidm)));

% ECC
[~,nt,nbCoils]=size(fidm);
fid = complex(zeros(size(fidm)));
for ical=1:nbCoils
    for jcal=1:nt
        if nbCoils==1
            % 2D matrix
            tempM =  fidm(:,jcal,ical);
            tempW = fidw(:,jcal);
            phiM = angle(tempM) - angle(tempM(1));
            phiW = angle(tempW) - angle(tempW(1));
            fid(:,jcal,ical) = abs(tempM).*(exp(1i*(phiM-phiW)));
            %fid(:,jcal,ical)=fidm(:,jcal,ical)./exp(1i*angle(fidw(:,jcal)));
        else
            % 3D matrix
            tempM =  fidm(:,jcal,ical);
            tempW = fidw(:,ical);
            phiM = angle(tempM) - angle(tempM(1));
            phiW = angle(tempW) - angle(tempW(1));
            fid(:,jcal,ical) = abs(tempM).*(exp(1i*(phiM-phiW)));
            %fid(:,jcal,ical)=fidm(:,jcal,ical)./exp(1i*angle(fidw(:,ical)));
        end
    end
end

% Redo DC correction
fid=fid-mean(fid(length(fid)-200:length(fid)));

return
