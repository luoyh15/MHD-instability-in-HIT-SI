function p = p_psi(psi)
%% the parameters determine the equilibrium
global R E a b T0 Psi_s

p = -4*Psi_s*(1+E^2)*(psi-Psi_s)/(a^2*E^2*R^2);