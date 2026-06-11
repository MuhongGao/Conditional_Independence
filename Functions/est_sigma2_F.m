function [sigma2_est] = est_sigma2_F(Z, Y, q ,p)

% X_Z = [X,Z];

n = length(Y);

k_plus_1 = 4;
[idx, ~] = knnsearch(Z, Z, 'K', k_plus_1);
I_x_all = idx(:, 2:k_plus_1);     % 去掉自身
%dist = dist(:, 2:end);

I_x_1 = I_x_all(:,1);
I_x_2 = I_x_all(:,2);
I_x_3 = I_x_all(:,3);


% [idx_2, ~] = knnsearch(X_Z, X_Z, 'K', 2);
% I_xz_1 = idx_2(:, 2);
%[I_xz_1,dist_2] = k_nearest_neighbors(X_Z, 1);

%X_1 = X(I_x,:);
R_Y = ceil(tiedrank(Y));
%R_Y = maxrank(Y);

R_Ni_1 = R_Y(I_x_1);
%R_Mi_1 = R_Y(I_xz_1);

R_Ni_2 = R_Y(I_x_2);
R_Ni_3 = R_Y(I_x_3);


a = min(R_Y, R_Ni_1);
a2 = min(R_Y, R_Ni_2);
a3 = min(R_Ni_2, R_Ni_3);

%b = min(R_Y, R_Mi_1);

% zeta = sum(b - a)/n^2;
% chi = sum(R_Y - a)/n^2;
% 
% xi = zeta/chi;

%qm = [0.67, 0.62, 0.59, 0.57];
qm = [0.667, 0.622, 0.593,0.573,0.558,0.547,0.538,0.531,0.528,0.521];

%om = [0.5, 0.618, 0.71, 0.76];
om = [0.5, 0.633, 0.709, 0.763, 0.805,0.84, 0.871, 0.898,0.923,0.946];

var_est = sum(a.*(a-a3))/n^3;
cov_est = sum(a.*(a2-a3))/n^3;



sigma2_est = (2 + qm(q) + qm(q+p)) * var_est + ...
    (om(q) + om(q+p) - 2*qm(q) - 2*qm(q+p) - 4) * cov_est;




end