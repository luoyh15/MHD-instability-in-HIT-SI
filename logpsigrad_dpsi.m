function dlogpsigdpsi = logpsigrad_dpsi(r,z)
global psi_dr_rz psi_dz_rz  psi_grad_norm...
    psi_grad_square_dr_rz psi_grad_square_dz_rz

dlogpsigdpsi = 1./(2*psi_grad_norm(r,z).^3).*...
    (psi_grad_square_dr_rz(r,z).*psi_dr_rz(r,z)+...
    psi_grad_square_dz_rz(r,z).*psi_dz_rz(r,z));

