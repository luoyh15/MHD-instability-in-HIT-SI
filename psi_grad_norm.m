function gradpsi = psi_grad_norm(r,z)

% the derivatives of psi
global psi_dr_rz psi_dz_rz
dpsidr = psi_dr_rz(r,z);
dpsidz = psi_dz_rz(r,z);
% the norm of gradient of psi
gradpsi = sqrt(dpsidr.^2+dpsidz.^2);
