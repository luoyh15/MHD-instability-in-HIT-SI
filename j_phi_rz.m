function j = j_phi_rz(r,z)
% psi
psi = psi_rz(r,z);
%% toroidal current density
j = r.*p_dpsi(psi)+1./(2*r).*2.*T_psi(psi).*T_dpsi(psi);