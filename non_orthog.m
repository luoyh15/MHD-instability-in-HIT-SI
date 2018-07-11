function beta_chi = non_orthog(s,chi)
global Psi_s psi_grad_square_dr_rz psi_grad_square_dz_rz psi_dr_rz psi_dz_rz

% dimension 
[n_s, n_chi] = size(s);
% initialize
beta_chi = zeros([n_s,n_chi]);
for i = 1:n_s
    % psi
    psi_i = s(i,1)^2*Psi_s;
    % the integral path
    chi_path = linspace(0,pi,200);
    s_path = ones(size(chi_path))*s(i,1);
    % the corresponding (r,z)
    [r_path,z_path] = schi2rz(s_path,chi_path);
    % the gradient of psi
    psi_grad_path = psi_grad_norm(r_path,z_path);
    % the derivative of psi gradient square in the normal direction of
    % constant psi surface
    psi_grad2_n_path = (psi_grad_square_dr_rz(r_path,z_path)...
        .*psi_dr_rz(r_path,z_path)+psi_grad_square_dz_rz(r_path,z_path)...
        .*psi_dz_rz(r_path,z_path))./(psi_grad_path);
    % the length of dchi
    dchi_path = [diff(chi_path),0];
    % the kernel of integral over consant psi path
    kernal_path = (j_phi_rz(r_path,z_path)./psi_grad_path+...
        1./T_psi(psi_i).*T_dpsi(psi_i)-...
        psi_grad2_n_path./psi_grad_path.^2-...
        1./q_psi(psi_i).*q_dpsi(psi_i)).*dchi_path;
    % integral over the path to get non-orthogonality
    beta_path = zeros(size(chi_path));
    for j =1:length(chi_path)
        beta_path(j) = sum(kernal_path(1:j));
    end
    beta_chi(i,:) = interp1(chi_path,beta_path,chi(i,:));
end
    

