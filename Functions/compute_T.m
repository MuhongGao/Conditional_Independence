function [T, S1, S2, tau, kappa] = compute_T(Z, Y, X, mode)
%
%  Function: compute Tn, kappa_n, tau_n

n = length(Y);

Z_X = [Z,X];

[idx, ~] = knnsearch(Z, Z, 'K', 2);

I_z_1 = idx(:,2);

[idx_2, ~] = knnsearch(Z_X, Z_X, 'K', 2);
I_zx_1 = idx_2(:, 2);

R_Y = ceil(tiedrank(Y));

R_Ni_1 = R_Y(I_z_1);
R_Mi_1 = R_Y(I_zx_1);

a = min(R_Y, R_Ni_1);
b = min(R_Y, R_Mi_1);

tau = sum(b - a)/n^2;

% mode = 0: only compute tau;
% mode = 1: compute all output variables.

T = 0;
kappa = 0;
S1 = 0;
S2 = 0;

if mode == 1

    kappa = sum(R_Y - a)/n^2;

    T = tau/kappa;

    S1 = sum(b)/n^2;

    S2 = sum(a)/n^2;

end

end