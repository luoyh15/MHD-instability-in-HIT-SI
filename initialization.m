%% set global variable
global R E a b T0 Psi_s lamda n gamma r_min r_max z_min z_max ...
    psi_dr_rz psi_dz_rz q0  psi_grad_norm ...
    psi_grad_square_dr_rz psi_grad_square_dz_rz

%% the parameters determine the equilibrium
R = 1; E = 2; a = 1/3; b = 0; T0 = 1;
Psi_s = 1; %the psi value on the plasma surface
lamda = 1; n = 2; gamma = 5/3;
%% the boundary of plasma
r_min = R*sqrt(1-2*a); r_max = R*sqrt(1+2*a);
r_temp = R*(1-4*a^2)^(1/4);
z_min = -sqrt((4*R^4*a^2-(r_temp^2-R^2)^2)*E^2/(4*r_temp^2));
z_max = -z_min;
%% the derivatives of psi
syms rr zz
fpsi = @(rr,zz) psi_rz(rr,zz);
psi_dr_rz = matlabFunction(diff(fpsi,rr));
psi_dz_rz = matlabFunction(diff(fpsi,zz));
%% the norm of psi gradient
fpsigradnorm = @(rr,zz) sqrt(psi_dr_rz(rr,zz).^2+psi_dz_rz(rr,zz).^2);
psi_grad_norm = matlabFunction(fpsigradnorm(rr,zz));
%% the derivatives of the square of psi gradient
psi_grad_square_dr_rz = matlabFunction(diff(psi_grad_norm(rr,zz).^2,rr));
psi_grad_square_dz_rz = matlabFunction(diff(psi_grad_norm(rr,zz).^2,zz));
% %% derivative of T
% syms pp
% fT = @(pp) T_psi(pp);
% T_dpsi = matlabFunction(diff(fT,pp));
% %% derivative of P
% fp = @(pp) p_psi(pp);
% p_dpsi = matlabFunction(diff(fp,pp));
%% safety factor at magnetic axis
q0 = q_psi(psi_rz(R,0));
