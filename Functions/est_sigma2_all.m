function [sigma2_est] = est_sigma2_all(Z, X, Y, T)
%
% Function: estimate sigma^2 = lim var(\tT_n), under general dependence.
%
%----------------------------
% Values of qd, od 
qd = [0.667, 0.622, 0.593,0.573,0.558,0.547,0.538,0.531,0.528,0.521];


od = [0.5, 0.633, 0.709, 0.763, 0.805,0.84, 0.871, 0.898,0.923,0.946];
%----------------------------


[n,q] = size(Z);
% 
[~,p] = size(X);

W = [X,Z];

k_plus_1 = 4;
[idz, ~] = knnsearch(Z, Z, 'K', k_plus_1);
I_z_all = idz(:, 2:k_plus_1);     % exclude self


I_z_1 = I_z_all(:,1);
I_z_2 = I_z_all(:,2);
I_z_3 = I_z_all(:,3);


[idw, ~] = knnsearch(W, W, 'K', k_plus_1);
I_w_all = idw(:, 2:k_plus_1);     % exclude self

I_w_1 = I_w_all(:,1);
I_w_2 = I_w_all(:,2);
I_w_3 = I_w_all(:,3);


%=============================================

R_Y = ceil(tiedrank(Y));
%R_Y = maxrank(Y);

R_Ni_1 = R_Y(I_z_1);
R_Ni_2 = R_Y(I_z_2);
R_Ni_3 = R_Y(I_z_3);

R_Mi_1 = R_Y(I_w_1);
R_Mi_2 = R_Y(I_w_2);
R_Mi_3 = R_Y(I_w_3);


a = min(R_Y, R_Ni_1);
a2 = min(R_Y, R_Ni_2);
a3 = min(R_Ni_2, R_Ni_3);
a4 = min(R_Ni_1, R_Ni_2);


b = min(R_Y, R_Mi_1);
b2 = min(R_Y, R_Mi_2);
b3 = min(R_Mi_2, R_Mi_3);
b4 = min(R_Mi_1, R_Mi_2);


a_sort = sort(a);

R_Ya = discretize(R_Y-0.001, [-Inf; a_sort(:); Inf]);
% R_Ya : n by 1.  R1(i) = #{j: a(j) < R_Y(i)}, i.e., the rank of R_Y(i) among
% a(1),...  a(n). 
% with -0.001: <;  else, <=
count_Ya = n-R_Ya +1;
% count_Ya: n by 1, count_Ya(i) =  #{j: a(j) >= R_Y(i)}


b_sort = sort(b);

R_Yb = discretize(R_Y-0.001, [-Inf; b_sort(:); Inf]);
count_Yb = n - R_Yb +1;


R_ab = discretize(a, [-Inf; b_sort(:); Inf]);
% R_ab : n by 1.  R2(i) = #{j: b(j) <= a(i)}, i.e., the rank of a(i) among
% b(1),...  b(n). 

count_ab = n - R_ab +1;
% count_ab: n by 1, count_ab(i) =  #{j: b(j)> a(i)}

R_ba = discretize(b, [-Inf; a_sort(:); Inf]);
count_ba = n - R_ba +1;

%=======================================================
% 1. compute covariance: sigma_1,2
%=======================================================


U1 = sum(a.*b)/n^3;
U2 = sum(a.*b4)/n^3;
U3 = sum(b.*a4)/n^3;

U4 =  sum(b.*count_Ya)/n^3;
U5 =  sum(b4.*count_Ya)/n^3;

U6 = sum(a.*count_Yb)/n^3;
U7 = sum(a4.*count_Yb)/n^3;

U8 = (sum(count_ab.*a) + sum(count_ba.*b))/n^3;
U9 = (sum(a)/n^2)*(sum(b)/n^2);

U_all = [U1, U2, U3, U4, U5, U6, U7, U8, U9];

coef = [4, -2, -1, 2, -1, 2, -1, 1, -4];

cov_est = sum(U_all.*coef);


%=======================================================
% 2. compute sigma_2^2
%=======================================================



R_aa = discretize(a, [-Inf; a_sort(:); Inf]);
% R2 : n by 1.  R2(i) = #{j: a(j) <= a(i)}, i.e., the rank of a(i) among
% a(1),...  a(n). 

count_aa = n - R_aa +1;
% count_aa: n by 1, count_aa(i) =  #{j: a(j) > a(i)}


T1 = sum(a.*a)/n^3;
T2 = sum(a.*a2)/n^3;
T3 = sum(a.*a3)/n^3;

T4 = sum(a.*count_Ya)/n^3;
T5 = sum(a4.*count_Ya)/n^3;
T6 = 2*sum(count_aa.*a)/n^3;

T7 = (sum(a)/n^2)^2;

T_all = [T1,T2,T3,T4,T5,T6,T7];
d = q;
coef = [1+qd(d), 2-2*qd(d)+od(d), -(2-qd(d)+od(d)), 4, -2, 1, -4];

%================  used for refinement of T6 ========
% % % % [u,~,idx] = unique(a);
% % % % cnt = accumarray(idx,1);
% % % % S = sum(a.*cnt(idx));
% % % % 
% % % % T6 =  (2*sum(count_aa.*a)+S)/n^3;
% T6 in this part is identical to that in paper, but only with o(1)
% difference, thus, for simplicity, we don't include this refinement.
%============================================================


sigma2_2 = sum(T_all.*coef);


%=======================================================
% 3. compute sigma_1^2
%=======================================================


R_bb = discretize(b, [-Inf; b_sort(:); Inf]);

count_bb = n - R_bb +1;


S1 = sum(b.*b)/n^3;
S2 = sum(b.*b2)/n^3;
S3 = sum(b.*b3)/n^3;

S4 = sum(b.*count_Yb)/n^3;
S5 = sum(b4.*count_Yb)/n^3;
S6 = 2*sum(count_bb.*b)/n^3;
S7 = (sum(b)/n^2)^2;

S_all = [S1,S2,S3,S4,S5,S6,S7];
d = p+q;
coef2 = [1+qd(d), 2-2*qd(d)+od(d), -(2-qd(d)+od(d)), 4, -2, 1, -4];



sigma2_1 = sum(S_all.*coef2);


% Output result ---------------
sigma2_est = sigma2_1 + (1-T)^2 * sigma2_2 - 2*(1-T)* cov_est;

sigma2_est = max(0.001, sigma2_est);




end