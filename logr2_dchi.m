function dlogr2dchi = logr2_dchi(s,chi)
n_s = length(s);
n_chi = length(chi);
% the dchi
dchi = 1/(20*n_chi);
% central difference
[r1,~] = schi2rz(s,chi(2:end-1)-dchi);
[r2,~] = schi2rz(s,chi(2:end-1)+dchi);
dlogr2dchi = zeros([n_s,n_chi]);
dlogr2dchi(:,2:end-1) = (log(r2.^2)-log(r1.^2))/(2*dchi);
end