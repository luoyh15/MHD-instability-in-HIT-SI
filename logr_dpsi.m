function dlogrdpsi = logr_dpsi(r,z)
global psi_dr_rz psi_dz_rz  psi_grad_norm
dlogrdpsi = 1./r.*psi_dr_rz(r,z)./psi_grad_norm(r,z);