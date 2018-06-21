function plot_equilibrium()
X = linspace(0.5,1.5,100);
Z = linspace(-0.5,0.5,100);
[X,Z] = meshgrid(X,Z);
psi = cal_equilibrium(X,Z);
psiB = 
contour(X,Z,psi);