function psi = psi_rz(r,z)

%% the parameters determine the equilibrium
global R E a b Psi_s
%% magnetic flux profile or the psi coordinate
psi = Psi_s/(a^2*R^4)*((r.^2+b^2).*z.^2/E^2+(r.^2-R^2).^2/4);
