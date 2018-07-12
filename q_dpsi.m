function dqdpsi = q_dpsi(psi)
global Psi_s
% psi for fit function
psi_step = Psi_s/1000;
n_fit = max(10,length(psi)*5);
psi_min = max(0,min(psi)/2);
psi_max = min(max(psi)*2,Psi_s);
psi_fit = linspace(psi_min,psi_max,n_fit);
q_fit = q_psi(psi_fit);
dq_fit = [0,diff(q_fit)./diff(psi_fit)];
dq_fit(1) = dq_fit(2);

dqdpsi = interp1(psi_fit,dq_fit,psi);
end