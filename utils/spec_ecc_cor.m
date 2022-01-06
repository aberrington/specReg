function [outfid]= spec_ecc_cor(infid,inref,lsfid_m,lsfid_w)
         phi_w_ref = 0.0;
         phi_m_ref = 0.0;
        
         np_metab=size(infid,1);
         np_water=size(inref,1);
         array= size(infid,2);
         array_w=size(inref,2);
         
         %fprintf (1,'\n******* EDDY CURRENT CORRECTION *******\n');
         %fprintf (1, 'np_metab = %d   np_water = %d\n',np_metab,np_water);
         %fprintf (1, 'lsfid(M) = %d   lsfid(W) = %d\n',lsfid_m,lsfid_w);
         if (lsfid_m > 0)
             lsfid_m=lsfid_m+1;
         infid_new=[infid(lsfid_m:end,:);zeros(lsfid_m-1,array)];
         else
             infid_new=infid;
         end
             
         if (lsfid_w > 0)
             lsfid_w=lsfid_w+1;
             inref_new=[inref(lsfid_w:end,:);zeros(lsfid_w-1,array_w)];
         else
             inref_new=inref;
         end
             
         
 
     %phi_w_ref=angle(inref_new(1,1));
     phi_w=angle(inref_new(:,1))-ones(np_water,1)*phi_w_ref;
     for i=1:array
         %phi_m_ref=angle(infid_new(1,i));
         phi_m_ref=0;
         phi_m=angle(infid_new(:,i))-ones(np_metab,1)*phi_m_ref-phi_w;
         abs_data = abs(infid_new(:,i));
         outfid (:,i)= complex(abs_data .* cos(phi_m), abs_data .* sin(phi_m));
     end
end
 