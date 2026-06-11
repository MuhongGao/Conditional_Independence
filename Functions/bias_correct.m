function [Lhat] = bias_correct(X,Y, lambda)
%
% Function: compute bias estimator Lhat for bias correction
%
%------------------------------------------------

[n,~] = size(X);

[idx, ~] = knnsearch(X, X, 'K', 2);

nn_idx = idx(:,2);

BX = quadratic_design_matrix(X); % n x d(d+3)/2

% i th row of Z: P1(Xi), P2(Xi).... PK(Xi), where K is the number of basis.
% Here, K = d(d+3)/2;

intercept = ones(n,1);

BX = [intercept, BX]; % n x (d(d+3)/2 +1 )


[~,b] = size(BX);   % b = d(d+3)/2 +1 

BtB = BX'*BX + n*lambda * eye(b);    % b * b
  
invBtB = BtB \ eye(b);    % b * b, equivalent to inv(BtB)


total = 0;

for j = 1: n
   t = Y(j);
   Ind_Yt = (Y>=t) ;
   beta_hat = invBtB * (BX' * Ind_Yt);
   
   G = BX * beta_hat;  % n*1 vector: (Ghat(X1,t), Ghat(X2,t)...Ghat(X_n,t))
   % where t = Y(j)
   % That is, corresponds to the fitted G at points X1,... X_n
   v = G.* G(nn_idx) - G.^2;
   total = total + sum(v) - v(j);
end


Lhat = total/(n *(n-1));

end




function Z = quadratic_design_matrix(X)
% QUADRATIC_DESIGN_MATRIX
% Build a quadratic (order-2) polynomial design matrix without intercept
%
% Input:
%   X : n x d matrix
%
% Output:
%   Z : n x p matrix
%       Columns correspond to:
%       1) linear terms:        X_j
%       2) interaction terms:   X_j * X_k, j < k
%       3) squared terms:       X_j^2

    [n, d] = size(X);                 % n: sample size, d: dimension

    %--------------------------------------------------
    % 1) Linear terms
    %--------------------------------------------------
    Z_linear = X;                     % n x d

    %--------------------------------------------------
    % 2) Pairwise interaction terms
    %--------------------------------------------------
    num_inter = d * (d - 1) / 2;
    Z_inter  = zeros(n, num_inter);   % n x (d(d-1)/2)

    idx = 1;
    for j = 1:d
        for k = j+1:d
            Z_inter(:, idx) = X(:, j) .* X(:, k);
            idx = idx + 1;
        end
    end

    %--------------------------------------------------
    % 3) Squared terms
    %--------------------------------------------------
    Z_square = X.^2;                  % n x d

    %--------------------------------------------------
    % Final design matrix
    %--------------------------------------------------
    Z = [Z_linear, Z_inter, Z_square];  % n x (d(d+3)/2)
end