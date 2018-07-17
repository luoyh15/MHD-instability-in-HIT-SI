%% set global variable
% global R E a b T0 Psi_s lamda n gamma r_min r_max z_min z_max ...
%     psi_dr_rz psi_dz_rz q0  psi_grad_norm ...
%     psi_grad_square_dr_rz psi_grad_square_dz_rz
%% the parameters determine the equilibrium
R = 1; E = 2; a = 1/3; b = 0; T0 = 1;
Psi_s = 1; %the psi value on the plasma surface
lamda = 1; n = 2; gamma = 5/3;
%% the boundary of plasma
r_min = R*sqrt(1-2*a); r_max = R*sqrt(1+2*a);
r_temp = R*(1-4*a^2)^(1/4);
z_min = -sqrt((4*R^4*a^2-(r_temp^2-R^2)^2)*E^2/(4*r_temp^2));
z_max = -z_min;
rho_max = sqrt(max(R-r_min,r_max-R)^2+z_max^2);
%% r z coordinate system 
%the dimensions of matrixes to fit
n_r = 1000;
n_z = 500;
r = linspace(r_min,r_max,n_r);
z = linspace(0,z_max,n_z);
[r,z] = meshgrid(r,z);
%% psi matrix for fit
M_psi = Psi_s/(a^2*R^4)*((r.^2+b^2).*z.^2/E^2+(r.^2-R^2).^2/4);
% the derivatives of psi
[M_psiDr,M_psiDz] = gradient(M_psi,r,z);
% the norm of psi gradient
M_psiGradNorm = sqrt(M_psiDr.^2+M_psiDz.^2);
% the derivatives of the square of psi gradient
[M_psiGrad2Dr,M_psiGrad2Dz] = gradient(M_psiGradNorm.^2,r,z);
% the gradient psi direction's derivative of psi gradient square
M_psiGrad2Dpsi = (M_psiGrad2Dr.*M_psiDr+M_psiGrad2Dz.*M_psiDz)./M_psiGradNorm;

%% psi as an independent variable for flux surface quantities
n_psi = 100;
psi = linspace(0,1,n_psi)*Psi_s;
%% T 
V_T = sqrt(T0^2-4*Psi_s*b^2./(a^2*R^2*E^2).*psi);
% derivative of T
V_TDpsi = gradient(V_T,psi);

%% pressure
V_p = -4*Psi_s*(1+E^2)*(psi-Psi_s)/(a^2*E^2*R^2);
% derivative of p
V_pDpsi = gradient(V_p,psi);

%% toroidal current density
M_jphi = r.*interp1(psi,V_pDpsi,M_psi)+...
    1./(2*r)*2.*interp1(psi,V_T,M_psi).*interp1(psi,V_TDpsi,M_psi);

%% safety factor
% integrate over the constant psi path
V_q = zeros(size(psi));
for i_psi = 2:length(psi)% exclude the magnetic axis situation
    psi_i = psi(i_psi);
    
    % the start and end point of the loop integration
    path_left = fzero(@(r) psi_rz(r,0)-psi_i,[r_min*0.9,R]);
    path_right = fzero(@(r) psi_rz(r,0)-psi_i,[R,r_max*1.1]);
    % discretization of theta in [0, pi]
    n_theta = 360;
    theta_path = linspace(0,pi,n_theta);   
    % the initialize vector of rho
    rho_path = zeros(size(theta_path));
    for i = 2:1:(n_theta-1)% exclude the start and end point at which z = 0
        rho_path(i) = fzero(@(rho) ...
            psi_rz(R+rho*cos(theta_path(i)),rho*sin(theta_path(i)))-psi_i,...
            [0,rho_max]);
    end
    % the corresponding r and z
    r_path = R+rho_path.*cos(theta_path);
    z_path = rho_path.*sin(theta_path);
    % the norm of gradient of psi
    psi_grad_norm_path = interp2(r,z,M_psiGradNorm,r_path,z_path);
    % length of path element
    dl_path = rho.*pi/(n_theta-1);
    % calculate q from loop integration
    V_q(i_psi) = 1/pi*sum(T_psi(psi_i)./(r_path.*psi_grad_norm_path).*dl_path);
end
q_fit =  polyfit(psi(2:end),V_q(2:end),10);
V_q(1) = polyval(q_fit,0);
% safety factor at magnetic axis
q0 = V_q(1);
% derivative of q
V_qDpsi = gradient(V_q,psi);

%% s chi coordinate system
% the dimensions 
n_s = 100;
n_chi = 100;
s = linspace(0,1,n_s);
chi = linspace(0,pi,n_chi);
[s,chi] = meshgrid(s,chi);
%% the corresponding r z coordinate
% initialize the r z
M_r = zeros(size(s));
M_z = zeros(size(s));
% integral through psi = const line to get chi at one point
for i = 1:n_s % the ith column of matrix s is constant 
    if(s(:,i) == 0)
        M_r(:,i) = R;
        M_z(:,i) = 0;
        continue;
    end
    % psi
    psi_i = s(1,i)^2*Psi_s;
    % the start and end point of integral path
    path_left = fzero(@(r) psi_rz(r,0)-psi_i,[r_min*0.9,R]);
    path_right = fzero(@(r) psi_rz(r,0)-psi_i,[R,r_max*1.1]);
    % discretization of theta in [0, pi]
    n_theta = 360;
    theta_path = linspace(0,pi,n_theta);   
    % the initialize vector of rho
    rho_path = zeros(size(theta_path));
    rho(1) = path_left;
    rho(end) = path_right;
    for i_path = 2:1:(n_theta-1)% exclude the start and end point at which z = 0
        rho_path(i) = fzero(@(rho) ...
            psi_rz(R+rho*cos(theta_path(i_path)),...
            rho*sin(theta_path(i_path)))-psi_i,...
            [0,rho_max]);
    end
    % the corresponding r and z
    r_path = R+rho_path.*cos(theta_path);
    z_path = rho_path.*sin(theta_path);
    % T
    T_path = interp1(psi,V_T,psi_i);
    % q 
    q_path = interp1(psi,V_q,psi_i);
    % the norm of gradient of psi
    psi_grad_norm_path = interp2(r,z,M_psiGradNorm,r_path,z_path);
    % length of path element
    dl_path = rho.*pi/(n_theta-1);
    % calculate integral kernel 
    kernel_path = T_path./(q_path.*r_path.*psi_grad_norm_path).*dl_path;
    % initialize of chi_path
    chi_path = zeros(size(theta_path));
    % calculate the chi value of this point
    for k = 2:(n_path-1)
        chi_path(k) = sum(kernel_path(1:k));
    end
    chi_path(end) = pi;
    M_r(:,i) = interp1(chi_path,r_path,chi(:,i));
    M_z(:,i) = interp1(chi_path,z_path,chi(:,i));  
end

%% non-orthgonolity
% initialize
M_betachi = zeros(size(s));
for i = 1:n_s
    % psi
    psi_i = s(1,i)^2*Psi_s;
    % discretization of chi in [0, pi]
    n_path = 5*n_chi;
    chi_path = linspace(0,pi,n_path);
    s_path = s(1,i)*ones(size(chi_path));
    % the corresponding (r,z)
    r_path = interp2(s,chi,M_r,s_path,chi_path);
    z_path = interp2(s,chi,M_z,s_path,chi_path);
    % T
    T_path = interp1(psi,V_T,psi_i);
    % psi derivative of T
    T_dpsi_path = interp1(psi,V_TDpsi,psi_i);
    % q 
    q_path = interp1(psi,V_q,psi_i);
    % psi derivative of q
    q_dpsi_path = interp1(psi,V_qDpsi,psi_i);
    % the gradient of psi
    psi_grad_norm_path = interp2(r,z,M_psiGradNorm,r_path,z_path);
    % the derivative of psi gradient square in gradient psi direction
    psi_grad2_dpsi_path = interp2(r,z,M_psiGrad2Dpsi,r_path,z_path);
    % the toroidal current density
    j_phi_path = interp2(r,z,M_jphi,r_path,z_path);
    % length of path element
    dchi_path = pi/(n_path-1);
    % the kernel of integral over consant psi path
    kernal = (j_phi_path./psi_grad_norm_path+...
        1./T_path.*T_dpsi_path-...
        psi_grad2_dpsi_path./psi_grad_norm_path.^2-...
        1./q_path.*q_dpsi_path).*dchi_path;
    % integral over the path to get non-orthogonality
    beta_path = zeros(size(chi_path));
    for j =1:length(chi_path)
        beta_path(j) = sum(kernal(1:j));
    end
    M_betachi(:,i) = interp1(chi_path,beta_path,chi(:,i));
end
