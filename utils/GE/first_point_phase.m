function [phase_ref,spec_out]=first_point_phase(spec_in)

fid_in=ifft(ifftshift(spec_in));

for t=1:size(spec_in,2)
    
    phase_ref(t)=angle(fid_in(1,t));
        
    Z=fid_in(:,t);
    
    phase_fid=angle(Z);
    phase_corr=phase_fid-phase_ref(t);
    R=abs(Z);
    fid_out(:,t) = R.*exp(1i*phase_corr);
       
    
end

spec_out=fftshift(fft(fid_out));









end