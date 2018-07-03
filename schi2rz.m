function [r_path,z_path] = schi2rz(s, chi)
global Psi_s R r_max z_max psi_dr_rz psi_dz_rz
if size(s) ~= size(chi)
    msg = 'error: the size of s and chi must be equal';
    error(msg);
end
% a constant used for calculate corresponding rho of a given theta
rho_max = sqrt(r_max^2+z_max^2);

% initialize the r z
n_s = size(s,1);
n_c = size(c,1);
r = zeros(size(n_s,n_c));
z = zeros(size(n_s,n_c));
%% integral through psi = const line to get chi at one point
for i = 1:length(s(:))
    if(s(i) == 0)
        r(i,:) = R;
        z(i,:) = 0;
        continue;
    end
    % psi
    psi_i = s(i)^2*Psi_s;
    % the start and end point of integral path
    path_left = fzero(@(r) psi_rz(r,0)-psi_i,[r_min*0.9,R]);
    path_right = fzero(@(r) psi_rz(r,0)-psi_i,[R,r_max*1.1]);
    % the integral over r = path_right to r = path_left
    n_path = n_c*20;
    r_path = linspace(path_left,path_right,n_path);
    % initialize of z_path chi_path
    z_path = zeros(size(r_path));
    chi_path = zeros(size(r_path));
    for j = 1:n_path
        % calculate the chi value of this point
        if(j == 1)
            chi_path(j) = 0;
        elseif(j == n_path)
            chi_path(j) = pi;
        else
            z_path(j) = fzero(@(z) psi_rz(r_path(j),z)-psi_i,[0,z_max*1.1]);
            
        end
    end
    % the derivatives of psi
    dpsidr = psi_dr_rz(r_path,z_path);
    dpsidz = psi_dz_rz(r_path,z_path);
    % the norm of gradient of psi
    psi_grad_path = sqrt(dpsidr.^2+dpsidz.^2);
    % length of path element
    %dl_path = sqrt(1+(dpsidr./dpsidz).^2).*(path_right-path_left)./(n_r-1);
    dl_path = [sqrt(diff(r_path).^2+diff(z_path).^2),0];
    % calculate integral kernel 
    kernel_path = T_psi(psi_i)./(q_psi(psi_i).*r_path.*psi_grad_path).*dl_path;
    for k = n_path
        chi_path(k) = sum(kernel_path(1:k));
    end
    
end
