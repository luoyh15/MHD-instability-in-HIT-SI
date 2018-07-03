function T = T_psi(psi)
%% the parameters determine the equilibrium
global R E a b T0 Psi_s

T = sqrt(T0^2-4*Psi_s*b^2./(a^2*R^2*E^2).*psi);