function dTqds = T_q_ds(s)
global Psi_s
n_s = length(s);
ds = 1/(20*n_s);
% central difference
psi1 = (s-ds).^2*Psi_s;
T1 = T_psi(psi1);
q1 = q_psi(psi1);
psi2 = (s+ds).^2*Psi_s;
T2 = T_psi(psi2);
q2 = q_psi(psi2);
dTqds = (T2./q2-T1./q1)/(2*ds);
end
