function chi = chi_rz_constpsi(r,psi)
global R r_max r_min z_max psi_dr_rz psi_dz_rz

%% integral through psi = const line to get chi

% the integral from r = r_left to r_right to get all the chi in the line
% psi = const
% the start and end point of the integration
r_left = fzero(@(r) psi_rz(r,0)-psi,[r_min*0.9,R]);
r_right = fzero(@(r) psi_rz(r,0)-psi,[R,r_max*1.1]);
n_r = max(500,length(r)*5);
r_path = linspace(r_right,r_left,n_r);
% the initialize vector of z corrensponding to r
z_path = zeros(size(r_path));
for i_r = 2:1:(n_r-1)% exclude the start and end point at which z = 0
    z_path(i_r) = fzero(@(z) psi_rz(r_path(i_r),z)-psi,[0,z_max*1.1]);
end

dl_path = [sqrt(diff(r_path).^2+diff(z_path).^2),0];
dl_path(end) = dl_path(end-1);

dpsidr_path = psi_dz_rz(r_path,z_path);
dpsidz_path = psi_dr_rz(r_path,z_path);
psi_grad_path = sqrt(dpsidr_path.^2+dpsidz_path.^2);
chi_kernel = T_psi(psi)./(q_psi(psi).*r_path.*psi_grad_path).*dl_path;
% chi value on every point on the path
chi_path = zeros(size(chi_kernel));
for i = 2:length(chi_kernel)-1
    chi_path(i) = sum(chi_kernel(1:i));
end
chi_path(end) = pi;
chi = interp1(r_path,chi_path,r);

