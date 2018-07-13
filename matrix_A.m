function A = matrix_A(n_s,n_chi)
%% s and chi 
s = 0:1/n_s:1;
chi = -pi/(2*(n_chi-1)):pi/(n_chi-1):pi+pi/(2*(n_chi-1));
%% initialize matrix A
n_x = (3*n_s+1)*(2*n_chi+2);% the totle number of unknow variables
A = zeros(n_x,n_x);
%% matrix form of the seven dependent components
%initialize
XR = zeros(n_s*n_chi,n_x);
XI = zeros(n_s*n_chi,n_x);
DXDchiR = zeros(n_s*n_chi,n_x);
DXDchiI = zeros(n_s*n_chi,n_x);
DXDsR = zeros(n_s*n_chi,n_x);
DXDsI = zeros(n_s*n_chi,n_x);
YR = zeros(n_s*n_chi,n_x);
YI = zeros(n_s*n_chi,n_x);
DYDchiR = zeros(n_s*n_chi,n_x);
DYDchiI = zeros(n_s*n_chi,n_x);
VR = zeros(n_s*n_chi,n_x);
VI = zeros(n_s*n_chi,n_x);
DVDchiR = zeros(n_s*n_chi,n_x);
DVDchiI = zeros(n_s*n_chi,n_x);
% the real part of the matrixes
for i = 1:n_s
    for j = 1:n_chi
        % the coefficients of the discretization form
        XR_coef = [1, 1, 1, 1]/4;
        XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*(j-1)+1) = XR_coef(1);
        XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*j+1) = XR_coef(2);
        XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*(j-1)+1) = XR_coef(3);
        XR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*j+1) = XR_coef(4);
       
        DXDchiR_coef = [-1, -1, 1, 1]/(2*(chi(j+1)-chi(j)));
        DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*(j-1)+1) = DXDchiR_coef(1);
        DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*j+1) = DXDchiR_coef(2);
        DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*(j-1)+1) = DXDchiR_coef(3);
        DXDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*j+1) = DXDchiR_coef(4);
        
        DXDsR_coef = [-1, -1, 1, 1]/(2*(s(i+1)-s(i)));
        DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*(j-1)+1) = DXDsR_coef(1);
        DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-3)+2*j+1) = DXDsR_coef(2);
        DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*(j-1)+1) = DXDsR_coef(3);
        DXDsR((i-1)*n_chi+j,(2*n_chi+2)*(3*i)+2*j+1) = DXDsR_coef(4);
        
        YR_coef = [1, 1]/2;
        YR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*(j-1)+1) = YR_coef(1);
        YR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*j+1) = YR_coef(2);
        
        DYDchiR_coef = [-1, 1]/(chi(j+1)-chi(j));
        DYDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*(j-1)+1) = DYDchiR_coef(1);
        DYDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-2)+2*j+1) = DYDchiR_coef(2);
        
        VR_coef = [1, 1]/2;
        VR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*(j-1)+1) = VR_coef(1);
        VR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*j+1) = VR_coef(2);
        
        DVDchiR_coef = [-1, 1]/(chi(j+1)-chi(j));
        DVDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*(j-1)+1) = DVDchiR_coef(1);
        DVDchiR((i-1)*n_chi+j,(2*n_chi+2)*(3*i-1)+2*j+1) = DVDchiR_coef(2);
       
    end
end
% the image part of the seven dependent variables
XI(:,2:2:end) = XR(:,1:2:end);
DXDchiI(:,2:2:end) = DXDchiR(:,1:2:end);
DXDsI(:,2:2:end) = DXDsR(:,1:2:end);
YI(:,2:2:end) = YR(:,1:2:end);
DYDchiI(:,2:2:end) = DYDchiR(:,1:2:end);
VI(:,2:2:end) = VR(:,1:2:end);
DVDchiI(:,2:2:end) = DVDchiR(:,1:2:end);

%% the normalized equilibrium quantities in the middle of each cell
global Psi_s q0 gamma n
% s and chi the physical components
s = (1/(2*n_s):1/n_s:1-1/(2*n_s))';
chi = (0:pi/(n_chi-1):pi)';
% psi
psi = s.^2*Psi_s;
% r z
[r,z] = schi2rz(s, chi);
% the norm of gradient of psi
psi_grad = psi_grad_norm(r,z);
% mass density
rho = 1;
% plasma pressure 
p = gamma*p_psi(psi)/(q0*Psi_s);
% flux function
T = T_psi(psi);
% toroidal current density
jphi = j_phi_rz(r,z);
% safety factor
q = q_psi(psi);
% poloidal magnetic field
Bp = psi.*r.^2./(q0*psi_grad.^2);
% non-orthogonality
beta_chi = non_orthog(s,chi);
% dlog(r^2)/ds and dlog(r^2)/dchi
dlogr2ds = logr2_ds(s,chi);
dlogr2dchi = logr2_dchi(s,chi);
% H defined at the ERATO paper
H = 2*jphi.*psi.*r./(s.*psi_grad.^2)+q./T.*T_q_ds(s);
% K defined at the ERATO paper
K = 2*psi/q0.*(jphi.^2./psi_grad.^2+jphi./r.*logpsigrad_dpsi(r,z)...
    -p_dpsi(psi).*logr_dpsi(r,z));

% reshape every quantity to one column
s = reshape(repmat(s,1,5)',[],1);
r = reshape(r',[],1);
psi = reshape(repmat(psi,1,5)',[],1);
psi_grad = reshape(psi_grad,[],1);
p = reshape(repmat(p,1,5)',[],1);
T = reshape(repmat(T,1,5)',[],1);
q = reshape(repmat(q,1,5)',[],1);
beta_chi = reshape(beta_chi',[],1);
dlogr2ds = reshape(dlogr2ds',[],1);
dlogr2dchi = reshape(dlogr2dchi',[],1);
H = reshape(H',[],1);
K = reshape(K',[],1);

%% the matrix form of Is
I1R = 1./q.*DXDchiR-n*XI;
I1I = 1./q.*DXDchiI+n*XR;
I2R = DXDsR+DVDchiR;
I2I = DXDsI+DVDchiI;
I3R = H.*XR-beta_chi.*n.*q.*XI-n*q.*VI+beta_chi.*DXDchiR+DXDsR;
I3I = H.*XI+beta_chi.*n.*q.*XR+n*q.*VR+beta_chi.*DXDchiI+DXDsI;
I4R = dlogr2ds.*XR+dlogr2dchi.*(VR+YR)-n*q.*YI+DXDsR+DVDchiR+DYDchiR;
I4I = dlogr2ds.*XI+dlogr2dchi.*(VI+YI)+n*q.*YR+DXDsI+DVDchiI+DYDchiI;
I5R = XR;
I5I = XI;
%% the coefficient a b c d e f g h
J = q.*r.^2./T;
a = 2*q.^2.*psi.*r.^4./(J.^3.*psi_grad.^2);
b = T.^2.*r.^2./(2*Psi_s*J);
c = psi_grad.^2.*r.^2./(2*Psi_s*J);
d = r.^4*gamma.*p./(2*Psi_s*J);
e = K.*2.*r.^4*q0./J;
f = 2*rho.*psi.*T.*r.^2./(q.*psi_grad.^2);
g = rho.*psi_grad.^2.*q.*r.^4./(2*T*Psi_s);
h = rho.*r.^4.*T.*q/(2*Psi_s);
%% the matrix A
A = I1R'*(a./s.*I1R) + I1I'*(a./s.*I1I) + I2R'*(b./s.*I2R) ...
    + I2I'*(b./s.*I2I) + I3R'*(c./s.*I3R) + I3I'*(c./s.*I3I)...
    + I4R'*(d./s.*I4R) + I4I'*(d./s.*I4I)...
    - I5R'*(e./s.*I5R) - I5I'*(e./s.*I5I);
%% the matrix B
B = XR'*(f./s.*XR) + XI'*(f./s.*XI) +...
    (VR+YR-beta_chi.*XR)'*(g./s.*(VR+YR-beta_chi.*XR))+...
    (VI+YI-beta_chi.*XI)'*(g./s.*(VI+YI+beta_chi.*XI))+...
    YR'*(h./s.*YR) + YI'*(h./s.*YI);

%% symmetry conditions
U = diag(ones(n_x,1));
for i = 1:n_s+1
    % symmetry conditions of XR
    U((2*n_chi+2)*(3*i-2)+1,(2*n_chi+2)*(3*i-2)+1) = 1; 
    U((2*n_chi+2)*(3*i-2)+1,(2*n_chi+2)*(3*i-2)+3) = 1;
    U((2*n_chi+2)*(3*i-2)+2*n_chi+1,(2*n_chi+2)*(3*i-2)+2*n_chi-1) = 1; 
    U((2*n_chi+2)*(3*i-2)+2*n_chi+1,(2*n_chi+2)*(3*i-2)+2*n_chi+1) = 1;
    % symmetry conditions of XI
    U((2*n_chi+2)*(3*i-2)+2,(2*n_chi+2)*(3*i-2)+2) = 1; 
    U((2*n_chi+2)*(3*i-2)+2,(2*n_chi+2)*(3*i-2)+4) = -1;
    U((2*n_chi+2)*(3*i-2)+2*n_chi+2,(2*n_chi+2)*(3*i-2)+2*n_chi) = 1; 
    U((2*n_chi+2)*(3*i-2)+2*n_chi+2,(2*n_chi+2)*(3*i-2)+2*n_chi+2) = -1;
    % symmetry conditions of VR
    U((2*n_chi+2)*(3*i-1)+1,(2*n_chi+2)*(3*i-1)+1) = 1; 
    U((2*n_chi+2)*(3*i-1)+1,(2*n_chi+2)*(3*i-1)+3) = -1;
    U((2*n_chi+2)*(3*i-1)+2*n_chi+1,(2*n_chi+2)*(3*i-1)+2*n_chi-1) = 1; 
    U((2*n_chi+2)*(3*i-1)+2*n_chi+1,(2*n_chi+2)*(3*i-1)+2*n_chi+1) = -1;
    % symmetry conditions of VI
    U((2*n_chi+2)*(3*i-1)+2,(2*n_chi+2)*(3*i-1)+2) = 1; 
    U((2*n_chi+2)*(3*i-1)+2,(2*n_chi+2)*(3*i-1)+4) = 1;
    U((2*n_chi+2)*(3*i-1)+2*n_chi+2,(2*n_chi+2)*(3*i-1)+2*n_chi) = 1; 
    U((2*n_chi+2)*(3*i-1)+2*n_chi+2,(2*n_chi+2)*(3*i-1)+2*n_chi+2) = 1;
    % symmetry conditions of YR
    U((2*n_chi+2)*(3*i)+1,(2*n_chi+2)*(3*i)+1) = 1; 
    U((2*n_chi+2)*(3*i)+1,(2*n_chi+2)*(3*i)+3) = -1;
    U((2*n_chi+2)*(3*i)+2*n_chi+1,(2*n_chi+2)*(3*i)+2*n_chi-1) = 1; 
    U((2*n_chi+2)*(3*i)+2*n_chi+1,(2*n_chi+2)*(3*i)+2*n_chi+1) = -1;
    % symmetry conditions of YI
    U((2*n_chi+2)*(3*i)+2,(2*n_chi+2)*(3*i)+2) = 1; 
    U((2*n_chi+2)*(3*i)+2,(2*n_chi+2)*(3*i)+4) = 1;
    U((2*n_chi+2)*(3*i)+2*n_chi+2,(2*n_chi+2)*(3*i)+2*n_chi) = 1; 
    U((2*n_chi+2)*(3*i)+2*n_chi+2,(2*n_chi+2)*(3*i)+2*n_chi+2) = 1;
end
        
    
    

