function [sigma2_est_b] = m_n_bootstrap(Z, Y, X, bootstrap_num, m)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here


n = length(Y);

tau_all = zeros(bootstrap_num, 1);

for i = 1: bootstrap_num

idx = randperm(n, m);
Zs = Z(idx,:);
Ys = Y(idx,:);
Xs = X(idx,:);

[~, ~, ~, tau, ~] = compute_T(Zs, Ys, Xs, 0);

%[zeta, ~] = compute_T(Zs, Ys, Zs);

tau_all(i) = tau;

end

sigma2_est_b = m * var(tau_all);


end