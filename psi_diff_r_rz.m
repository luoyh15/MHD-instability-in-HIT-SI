function psi_dr = psi_diff_r_rz(r,z)
syms rr zz
fpsi = @(rr,zz) psi_rz(rr,zz);
fpsi_dr = matFunction(diff(fpsi,rr));
psi_dr = fpsi_dr(r,z);
