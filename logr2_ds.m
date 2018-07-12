function dlogr2ds = logr2_ds(s,chi)
n_s = length(s);
n_chi = length(chi);
% the ds
ds = 1/(20*n_s);
% central difference
[r1,~] = schi2rz(s-ds,chi);
[r2,~] = schi2rz(s+ds,chi);
dlogr2ds = (log(r2.^2)-log(r1.^2))/(2*ds);
end
