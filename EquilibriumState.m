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
%% the dimensions of matrixes to fit
n_r = 100;
n_z = 100;
r = linspace(r_min,r_max,n_r);
z = linspace(0,z_max,n_z);
[r,z] = meshgrid(r,z);
%% psi matrix for fit
M_psi = Psi_s/(a^2*R^4)*((r.^2+b^2).*z.^2/E^2+(r.^2-R^2).^2/4);

%% the derivatives of psi
[M_psiDr,M_psiDz] = gradient(M_psi,r,z);

%% the norm of psi gradient
M_psiGradNorm = sqrt(M_psiDr.^2+M_psiDz.^2);

%% the derivatives of the square of psi gradient
[M_psiGrad2Dr,M_psiGrad2Dz] = gradient(M_psiGradNorm.^2,r,z);

%% derivative of T
% syms pp
% fT = @(pp) T_psi(pp);
% T_dpsi = matlabFunction(diff(fT,pp));
% %% derivative of P
% fp = @(pp) p_psi(pp);
% p_dpsi = matlabFunction(diff(fp,pp));
%% safety factor at magnetic axis
q0 = q_psi(psi_rz(R,0));
