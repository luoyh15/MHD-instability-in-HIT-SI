function beta_chi = non_orthogonality(s,chi)
% the step length for difference
d_step = 1e-8;
% the corresponding (r,z)
[r,z] = schi2rz(s,chi);
% the gradient of psi determine the direction of nu(the subscript of
% beta_chi)
global psi_dr_rz psi_dz_rz Psi_s
dpsidr = psi_dr_rz(r,z);
dpsidz = psi_dz_rz(r,z);
psi_grad = psi_grad_norm(r,z);
% initialize non-orthogonality
beta_chi = zeros(size(s));

for i = 1:size(s,1)
    
    if(s(i,1)<0.1)
        psi = s(i,1)^2*Psi_s;
        psi1 = (s(i,1)+d_step)^2*Psi_s;
        dpsi = psi1-psi;
      
        chi1 = chi(i,:);
        rtemp = r(i,:)+dpsi*dpsidr(i,:)./psi_grad(i,:).^2;
        ztemp = z(i,:)+dpsi*dpsidz(i,:)./psi_grad(i,:).^2;
        psitemp = psi_rz(rtemp,ztemp);
        stemp = sqrt(psitemp/Psi_s);
        chi2 = chi_rz_constpsi(rtemp,psi1);
        
        beta_chi(i,:) = (chi2-chi1)./(stemp-s(i,:));
    elseif(1-s(i,:)<d_step) 
        psi1 = (s(i,1)-d_step)^2*Psi_s;
        psi = s(i,1)^2*Psi_s;
        dpsi = psi1-psi;
        
        rtemp = r(i,:)+dpsi*dpsidr(i,:)./psi_grad(i,:).^2;
%         ztemp = z(i,:)-d_step*dpsidz(i,:)./psi_grad(i,:);
        chi1 = chi_rz_constpsi(rtemp,psi1);
        chi2 = chi(i,:);
        
        beta_chi(i,:) = (chi2-chi1)./d_step;
    else
        psi1 = (s(i,1)-d_step)^2*Psi_s;
        psi2 = (s(i,1)+d_step)^2*Psi_s;
        psi = s(i,1).^2*Psi_s;
        dpsi1 = psi1-psi;
        dpsi2 = psi2-psi;
        
        rtemp = r(i,:)+dpsi1*dpsidr(i,:)./psi_grad(i,:).^2;
        ztemp = z(i,:)+dpsi1*dpsidz(i,:)./psi_grad(i,:).^2;
        chi1 = chi_rz_constpsi(rtemp,psi1);
        psitemp = psi_rz(rtemp,ztemp);
        
        rtemp = r(i,:)+dpsi2*dpsidr(i,:)./psi_grad(i,:).^2;
        ztemp = z(i,:)+dpsi2*dpsidz(i,:)./psi_grad(i,:).^2;
        psitemp = psi_rz(rtemp,ztemp);
        chi2 = chi_rz_constpsi(rtemp,psi2);
        
        beta_chi(i,:) = (chi2-chi1)./(2*d_step);
    end
        
end
    
