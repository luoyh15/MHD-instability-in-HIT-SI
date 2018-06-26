function [r,z,psi,pressure] = cal_equilibrium()
 %calculate the equilibrium
%% the parameters determine the equilibrium
R = 1; E = 2; a = 1/3; b = 0; T0 = 1;
Psi_s = 1; %the psi value on the plasma surface 
 lamda = 1; n = 2;
%% the boundary of plasma

Xmin = R*sqrt(1-2*a); Xmax = R*sqrt(1+2*a);
Xtemp = R*(1-4*a^2)^(1/4);
Zmin = -sqrt((4*R^4*a^2-(Xtemp^2-R^2)^2)*E^2/(4*Xtemp^2));
Zmax = -Zmin;
%% the psi
nX = 100; nZ = 100;
r = linspace(0.9*Xmin, 1.05*Xmax, nX);
z = linspace(1.1*Zmin, 1.1*Zmax, nZ);
[r,z] = meshgrid(r,z);
%magnetic flux profile or the psi coordinate
psi = Psi_s/(a^2*R^4)*((r.^2+b^2).*z.^2/E^2+(r.^2-R^2).^2/4);
%pressure profile
pressure = -4*Psi_s*(1+E^2)/(a^2*E^2*R^2)*(psi-Psi_s);
%flux of poloidal current
flux_jp = sqrt(T0^2-4*Psi_s*b^2/(a^2*R^2*E^2)*psi);
%safety factor q 
q = 

z_temp = sqrt(



%% plot the equilibrium state
%v = linspace(0,psiB,8);
v = [0:0.5:Psi_s,Psi_s];
figure(1);
hold on;
[C,h] = contour(r,z,psi,v);
h.LineColor = 'k';
h.ShowText = 'on';
contourf(r,z,pressure,v);
axis square;
hold off;

