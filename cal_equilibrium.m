function [r,z,psi,pressure] = cal_equilibrium()
%calculate the equilibrium
%% the parameters determine the equilibrium
R = 1; E = 2; a = 1/3; b = 0; T0 = 1;
Psi_s = 1; %the psi value on the plasma surface
lamda = 1; n = 2;
%% the boundary of plasma

r_min = R*sqrt(1-2*a); r_max = R*sqrt(1+2*a);
r_temp = R*(1-4*a^2)^(1/4);
z_min = -sqrt((4*R^4*a^2-(r_temp^2-R^2)^2)*E^2/(4*r_temp^2));
z_max = -z_min;
%% the discretized equilibrium state
n_r = 100; n_z = 100;
r = linspace(0.9*r_min, 1.05*r_max, n_r);
z = linspace(1.1*z_min, 1.1*z_max, n_z);
[rgrid,zgrid] = meshgrid(r,z);

%magnetic flux profile or the psi coordinate
psi = Psi_s/(a^2*R^4)*((rgrid.^2+b^2).*zgrid.^2/E^2+(rgrid.^2-R^2).^2/4);

% s coordinate
n_s = 100;
s = linspace(0,1,n_s);
n_in = sqrt(n_s);
Psi_in = s(n_in)^2*Psi_s;
index_in = find(psi<Psi_in);
%psi_in = psi(psi<Psi_in);
psi_fit_in = fit([rgrid(index_in).^2,zgrid(index_in).^2], psi(index_in),'poly33');
%poloidal current flux
T = @(r,z) sqrt(T0^2-4*Psi_s*b^2./(a^2*R^2*E^2).*psi_fit_in(r,z));
psi_diff = @(r,z) differentiate(psi_fit_in,[r.^2,z.^2]).*[2*r,2*z];
norm_psi_gradient = @(r,z) norm(differentiate(psi_fit_in,[r,z]));

q = zeros(size(s));
for i_s = 2:n_s %this should be n_in, but in this case the fit function is perfect in the whole region
    psi_i = s(i_s)^2*Psi_s;
    path_left = fzero(@(r) psi_fit_in(r.^2,0)-psi_i,[r_min,R]);
    path_right = fzero(@(r) psi_fit_in(r.^2,0)-psi_i,[R,r_max]);
    path_function = @(r) fzero(@(z) psi_fit_in(r.^2,z.^2)-psi_i,[0,z_max*1.1]);
    %psi_diff = @(r,z) differentiate(psi_fit_in,[r,z]);
    %calculate q path integration
    n_r_int = 100;
    r_pace_int = (path_left - path_right)/n_r_int;
    r_int = linspace(path_left,path_right,100)';
    for r_i = linspace(path_left,path_right,n_r_int)
        
        z_i = path_function(r_i);
        psi_diff_i = psi_diff(r_i,z_i);
        dz_dr = psi_diff_i(1,1)./psi_diff_i(2,1);
        dl_i = sqrt(1+dz_dr.^2);
        norm_psi_gradient_i = sqrt(psi_diff_i(1,1).^2+psi_diff_i(2,1).^2);
        q_temp = T(r_i,z_i)./(r_i.*norm_psi_gradient_i).*dl_i.*r_pace_int;
        q(i_s) = q(i_s)+q_temp;
    end
%     psi_diff_r = @(r) psi_diff(r,path_function(r));
%     path_length_unit = @(r) 
    q = integral(@(r) T./(r.*norm_psi_gradient(r,path_function(r))).*path_length_unit(r),path_left,path_right);
    
    
end

%q near the magenetic axis
for i_s = 1:n_in
    psi_i = s(i_s)^2*Psi_s;
    
    
end







%pressure profile
pressure = -4*Psi_s*(1+E^2)./(a^2*E^2*R^2).*(psi-Psi_s);
%flux of poloidal current
T = sqrt(T0^2-4*Psi_s*b^2./(a^2*R^2*E^2).*psi);
%safety factor
q = cal_q_rz(psi,T,r,z);
%poloidal angle chi coordinate


%safety factor q
% q =

% z_temp = sqrt(



%% plot the equilibrium state
%v = linspace(0,psiB,8);
v = [0:0.5:Psi_s,Psi_s];
figure(1);
hold on;
[C,h] = contour(r,z,psi,v);
[C1,h1]=contour(r,z,psi,Psi_s);
h.LineColor = 'k';
h.ShowText = 'off';
%contourf(r,z,pressure,v);
axis square;
plot(C(1,:),C(2,:),'w.');
hold off;
end

function q = cal_q_rz(psi,T,r,z)
[n_r,n_z] = size(psi);
[rgrid,zgrid] = meshgrid(r,z);

%the norm of gradient psi
[psi_x,psi_y] = gradient(psi,r,z);
norm_psi_gradient = sqrt(psi_x.^2+psi_y.^2);

%calculate q on every point
q = zeros(size(psi));
for ir = 1:n_r
    for iz = 1:n_z
        if psi(ir,iz)>1
            q(ir,iz) = NaN;
            continue;
        end
        cdata = contourc(r,z,psi,[-1,psi(ir,iz)]);
        s = contourdata(cdata);
        r_on_path = s.xdata;
        z_on_path = s.ydata;
        %psi_on_path = interp2(r,z,psi,r_on_path,z_on_path);
        norm_psi_gradient_on_path = interp2(rgrid,zgrid,norm_psi_gradient,r_on_path,z_on_path);
        T_on_path = interp2(rgrid,zgrid,T,r_on_path,z_on_path);
        arc_length_on_path = sqrt(gradient(r_on_path).^2+gradient(z_on_path).^2);
        q(ir,iz) = 1/(2*pi)*sum(T_on_path./(r_on_path.*norm_psi_gradient_on_path).*arc_length_on_path);
    end
end
end





