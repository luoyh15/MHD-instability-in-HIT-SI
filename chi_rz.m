function chi = chi_rz(r,z)
global R r_max z_max psi_dr_rz psi_dz_rz
if size(r) ~= size(z)
    msg = 'error: the size of r and z must be equal';
    error(msg);
end
% a constant used for calculate corresponding rho of a given theta
rho_max = sqrt(r_max^2+z_max^2);
%% integral through psi = const line to get chi at one point
chi = zeros(size(r));
for i = 1:length(r(:))
    % psi
    psi_i = psi_rz(r(i),z(i));
    % polar coordinate
    rho_i = sqrt((r(i)-R).^2+z(i).^2);
    theta_i = acos((r(i)-R)./rho_i);
    if(theta_i==0||theta_i==pi)
        chi(i) = theta_i;
        continue;
    end  
    % the integral over theta = 0 to theta = theta_i
    n_theta = max(100,floor(theta_i/pi*500)+1);
    theta_path = linspace(0,theta_i,n_theta);
    rho_path = zeros(size(theta_path));
    for j = 1:length(theta_path)
        rho_path(j) = fzero(@(rho) psi_rz(R+rho.*cos(theta_path(j)),rho.*sin(theta_path(j)))-psi_i,...
        [0,rho_max]);
    end
    % the corresponding r, z, and dl
    r_path = R+rho_path.*cos(theta_path);
    z_path = rho_path.*sin(theta_path);
    dl_path = [sqrt(diff(r_path).^2+diff(z_path).^2),0];
    dl_path(end) = dl_path(end-1);
    % the rho derivative of chi
    psi_dr_path = psi_dz_rz(r_path,z_path);
    psi_dz_path = psi_dr_rz(r_path,z_path);
    psi_grad_path = sqrt(psi_dr_path.^2+psi_dz_path.^2);
    chi(i) = sum(T_psi(psi_i)./(q_psi(psi_i).*r_path.*psi_grad_path).*dl_path);
end

