function [X,Z,psi,P] = cal_equilibrium()
 %calculate the equilibrium
%% the parameters determine the equilibrium
R = 1; B0 = 1; E = 2; epsilon = 1/3;
q0 = 0.3; lamda = 1; n = 2;
%% the boundary of plasma
psiB = epsilon^2*pi*E*R^2*B0/q0;
Xmin = R*sqrt(1-2*epsilon); Xmax = R*sqrt(1+2*epsilon);
Xtemp = R*(1-4*epsilon^2)^(1/4);
Zmin = -sqrt((4*R^4*epsilon^2-(Xtemp^2-R^2)^2)*E^2/(4*Xtemp^2));
Zmax = -Zmin;
%% the psi
nX = 100; nZ = 100;
X = linspace(0.9*Xmin, 1.05*Xmax, nX);
Z = linspace(1.1*Zmin, 1.1*Zmax, nZ);
[X,Z] = meshgrid(X,Z);
psi = pi*B0/(E*R^2*q0)*(X.^2.*Z.^2+E^2/4*(X.^2-R^2).^2);

P = (1+E^2)*B0/(2*pi*E*R^2*q0)*(psiB-psi);
%v = linspace(0,psiB,8);
v = [0:0.5:psiB,psiB];
figure(1);
hold on;
[C,h] = contour(X,Z,psi,v);
h.LineColor = 'k';
h.ShowText = 'on';
contourf(X,Z,P,v);
axis square;
hold off;

