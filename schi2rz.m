function [r,z] = schi2rz(s, chi)
global Psi_s R r_min r_max z_max psi_dr_rz psi_dz_rz
if size(s) ~= size(chi)
    msg = 'error: the size of s and chi must be equal';
    error(msg);
end
% initialize the r z
n_s = size(s,1);
n_chi = size(s,2);
r = zeros(n_s,n_chi);
z = zeros(n_s,n_chi);
%% integral through psi = const line to get chi at one point
for i = 1:n_s
    if(s(i,1) == 0)
        r(i,:) = R;
        z(i,:) = 0;
        continue;
    end
    % psi
    psi_i = s(i,1)^2*Psi_s;
    % the start and end point of integral path
    path_left = fzero(@(r) psi_rz(r,0)-psi_i,[r_min*0.9,R]);
    path_right = fzero(@(r) psi_rz(r,0)-psi_i,[R,r_max*1.1]);
    % the integral over r = path_right to r = path_left
    n_path = max(500,10*n_chi);
    r_path = linspace(path_right,path_left,n_path);
    % initialize of z_path chi_path
    z_path = zeros(size(r_path));
    chi_path = zeros(size(r_path));
    % calculate the chi value of this point
    for j = 2:n_path-1
        z_path(j) = fzero(@(z) psi_rz(r_path(j),z)-psi_i,[0,z_max*1.1]);
    end
    % the norm of gradient of psi
    psi_grad_path = psi_grad_norm(r_path,z_path);
    % length of path element
    %dl_path = sqrt(1+(dpsidr./dpsidz).^2).*(path_right-path_left)./(n_r-1);
    dl_path = [sqrt(diff(r_path).^2+diff(z_path).^2),0];
    % calculate integral kernel 
    kernel_path = T_psi(psi_i)./(q_psi(psi_i).*r_path.*psi_grad_path).*dl_path;
    for k = 2:(n_path-1)
        chi_path(k) = sum(kernel_path(1:k));
    end
    chi_path(end) = pi;
    r(i,:) = interp1(chi_path,r_path,chi(i,:));
    z(i,:) = interp1(r_path,z_path,r(i,:));
    
end
