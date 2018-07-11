function A = matrix_A(s,chi)
%% s and chi should be same size (n_s,n_chi)  meshgrid type
[n_s,n_chi] = size(s);
%% initialize matrix A
n_x = (3*n_s+1)*(2*n_chi+2);% the totle number of unknow variables
A = zeros(n_x,n_x);
%% matrix form of the seven dependent components
%initialize
M_XR = zeros(n_s*n_chi,n_x);
M_XI = zeros(n_s*n_chi,n_x);
M_DXDchiR = zeros(n_s*n_chi,n_x);
M_DXDchiI = zeros(n_s*n_chi,n_x);
M_DXDsR = zeros(n_s*n_chi,n_x);
M_DXDsI = zeros(n_s*n_chi,n_x);
M_YR = zeros(n_s*n_chi,n_x);
M_YI = zeros(n_s*n_chi,n_x);
M_DYDchiR = zeros(n_s*n_chi,n_x);
M_DYDchiI = zeros(n_s*n_chi,n_x);
M_VR = zeros(n_s*n_chi,n_x);
M_VI = zeros(n_s*n_chi,n_x);
M_DVDchiR = zeros(n_s*n_chi,n_x);
M_DVDchiI = zeros(n_s*n_chi,n_x);
% the real part of the matrixes
for i = 1:n_s
    for j = 1:n_chi
        % the coefficients of the discretization form
        XR_coef = [1, 1, 1, 1]/4;
        M_XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*(j-1)+1) = XR_coef(1);
        M_XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*j+1) = XR_coef(2);
        M_XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*(j-1)+1) = XR_coef(3);
        M_XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*j+1) = XR_coef(4);
       
        DXDchiR_coef = [-1, -1, 1, 1]/(2*(chi(i,j+1)-chi(i,j)));
        M_DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*(j-1)+1) = DXDchiR_coef(1);
        M_DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*j+1) = DXDchiR_coef(2);
        M_DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*(j-1)+1) = DXDchiR_coef(3);
        M_DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*j+1) = DXDchiR_coef(4);
        
        DXDsR_coef = [-1, -1, 1, 1]/(2*(s(i+1,j)-s(i,j)));
        M_DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*(j-1)+1) = DXDsR_coef(1);
        M_DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*j+1) = DXDsR_coef(2);
        M_DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*(j-1)+1) = DXDsR_coef(3);
        M_DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*j+1) = DXDsR_coef(4);
        
        YR_coef = [1, 1]/4;
        M_YR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*(j-1)+1) = YR_coef(1);
        M_YR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*j+1) = YR_coef(2);
        
        DYDchiR_coef = [-1, 1]/(chi(i,j+1)-chi(i,j));
        M_DYDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*(j-1)+1) = DYDchiR_coef(1);
        M_DYDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*j+1) = DYDchiR_coef(2);
        
        VR_coef = [1, 1]/4;
        M_VR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*(j-1)+1) = VR_coef(1);
        M_VR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*j+1) = VR_coef(2);
        
        DVDchiR_coef = [-1, 1]/(chi(i,j+1)-chi(i,j));
        M_DVDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*(j-1)+1) = DVDchiR_coef(1);
        M_DVDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*j+1) = DVDchiR_coef(2);
       
    end
end
% the image part of the seven dependent variables
M_XI(:,2:2:end) = M_XR(:,1:2:end);
M_DXDchiI(:,2:2:end) = M_DXDchiR(:,1:2:end);
M_DXDsI(:,2:2:end) = M_DXDsR(:,1:2:end);
M_YI(:,2:2:end) = M_YR(:,1:2:end);
M_DYDchiI(:,2:2:end) = M_DYDchiR(:,1:2:end);
M_VI(:,2:2:end) = M_VR(:,1:2:end);
M_DVDchiI(:,2:2:end) = M_DVDchiR(:,1:2:end);

%% the information of equilibrium state 
global Psi_s q0 
psi = s.^2*Psi_s;
[r,z] = schi2rz(s, chi);
% the norm of gradient of psi
psi_grad = psi_grad_norm(r,z);
% safety factor
q = q_psi(psi);
% poloidal magnetic field
Bp = psi.*r.^2./(q0*psi_grad.^2);







