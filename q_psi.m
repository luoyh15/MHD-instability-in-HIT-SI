function q = q_psi(psi)
%% the parameters determine the equilibrium
global R E a b T0 lamda n r_min r_max z_min z_max psi_dr_rz psi_dz_rz
%% the loop integration
q = zeros(size(psi));
for i_psi = 1:length(psi(:))
    psi_i = psi(i_psi);
    
    % the magnetic axis situation
    if psi_i == psi_rz(R,0)
        psi_i = psi_i+1e-8;
    end
    
    % the start and end point of the loop integration
    path_left = fzero(@(r) psi_rz(r,0)-psi_i,[r_min*0.9,R]);
    path_right = fzero(@(r) psi_rz(r,0)-psi_i,[R,r_max*1.1]);
    % discretization of r in [path_left, path_right]
    n_r = 100;
    r = linspace(path_left,path_right,n_r);   
    % the initialize vector of z corrensponding to r
    z = zeros(size(r));   
    for i_r = 2:1:(n_r-1)% exclude the start and end point at which z = 0
        z(i_r) = fzero(@(z) psi_rz(r(i_r),z)-psi_i,[0,z_max*1.1]);
    end
    % the derivatives of psi
    dpsidr = psi_dr_rz(r,z);
    dpsidz = psi_dz_rz(r,z);
    % the norm of gradient of psi
    dpsi_norm = sqrt(dpsidr.^2+dpsidz.^2);
    % length of path element
    %dl_path = sqrt(1+(dpsidr./dpsidz).^2).*(path_right-path_left)./(n_r-1);
    dl_path = [sqrt(diff(r).^2+diff(z).^2),0];
    % calculate q from loop integration
    q(i_psi) = sum(1/(2*pi)*T_psi(psi_i)./(r.*dpsi_norm).*dl_path*2);
end
