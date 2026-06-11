% Name: Sim_ECP.m
% Function: Conduct simulations realted to ECPs in Section 5.1 (Tables 2 and 3)
%
% Input: Input the model number (1 or 2) and replication num, following the prompt.
%
% Note: This code performs parallel computing. Default Corenum = 64.
%       User may change the "Corenum" according to the machines capacity.
%
% Elasped time: about 50 min for rep_num = 1000. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('Functions/');

model = input('Input the model number (1 or 2): ');
rep_num = input('Input the replication number (100, 500, or 1000): ');

p_all = [1,1,3,3,5];
q_all = [1,2,1,3,5];
pq_num = length(p_all);

rho_all = [0,0.3,0.5, 0.7, 0.9];

rho_num = length(rho_all);

n_all = [1000,5000,10000];

n_num = length(n_all);

%rep_num = 1000;


T_all = zeros(pq_num, rho_num, n_num, rep_num);

Tbc_all = zeros(pq_num, rho_num, n_num, rep_num);

tau_all = zeros(pq_num, rho_num, n_num, rep_num);

kappa_all = zeros(pq_num, rho_num, n_num, rep_num);

sigma2_all = zeros(pq_num, rho_num, n_num, rep_num);


CI_all = zeros(2, pq_num, rho_num, n_num, rep_num);

CIbc_all = zeros(2, pq_num, rho_num, n_num, rep_num);



CoreNum = 64;


z_alpha_half = norminv(0.95);  %alpha = 0.1




if isempty(gcp('nocreate'))
    parpool(CoreNum);
end

tic;

parfor (i = 1: rep_num, CoreNum)

    if mod(i,100) == 1
        fprintf([' rep = ', num2str(i), '\n']);
    end

    p_all_c = p_all;
    q_all_c = q_all;
    rho_all_c = rho_all;
    n_allc = n_all;

    T_this = zeros(pq_num,rho_num, n_num);
    Tbc_this = zeros(pq_num,rho_num, n_num);
    kappa_this = zeros(pq_num,rho_num, n_num);
    tau_this = zeros(pq_num,rho_num, n_num);
    sigma2_this = zeros(pq_num,rho_num, n_num);
    CI_this = zeros(2, pq_num,rho_num, n_num);
    CIbc_this = zeros(2, pq_num,rho_num, n_num);


    for j = 1: pq_num

        p = p_all_c(j);
        q = q_all_c(j);

        for k = 1:rho_num

            rho = rho_all_c(k);

            if rho == 0
                H = 0;
            else
                H = 1;
            end

            for r = 1 : n_num

                n = n_allc(r);

                lambda = n^(-0.85);

                seed = i+1000;
                rng(seed, 'twister');


                % Generate data
                [Z,Y,X] = data_gener(n,q,p, model, H, rho);


                % Compute Tn, kappa_n, tau_n
                [T, S1, S2, tau, kappa] = compute_T(Z, Y, X, 1);


                % Bias correction -----------
                Lhat1 = bias_correct([Z,X],Y, lambda);

                Lhat2 = bias_correct(Z,Y, lambda);

                tau_bias_est = Lhat1 - Lhat2;

                kappa_bias_est = - Lhat2;

                tau_bc = tau - tau_bias_est;

                kappa_bc = kappa - kappa_bias_est;

                Tbc =  tau_bc/kappa_bc;
                % ------------------------------

                % Estimate limiting variance
                sigma2_est = est_sigma2_all(Z, X, Y, T);


                % Compute confidence intervales ------------
                CI_low = T - z_alpha_half / kappa *sqrt(sigma2_est/n);
                CI_up = T + z_alpha_half / kappa *sqrt(sigma2_est/n);

                CI = [CI_low, CI_up];

                CIbc_low = Tbc - z_alpha_half / kappa *sqrt(sigma2_est/n);
                CIbc_up = Tbc + z_alpha_half / kappa *sqrt(sigma2_est/n);

                CIbc = [CIbc_low, CIbc_up];
                % ------------------------------------------



                T_this(j,k,r) = T;
                Tbc_this(j,k,r) = Tbc;
                kappa_this(j,k,r) = kappa;
                tau_this(j,k,r) = tau;
                sigma2_this(j,k,r) = sigma2_est;

                CI_this(:,j,k,r) = CI;
                CIbc_this(:,j,k,r) = CIbc;
            end


        end
    end

    T_all(:,:,:,i) = T_this;
    Tbc_all(:,:,:,i) = Tbc_this;
    kappa_all(:,:,:,i) = kappa_this;
    tau_all(:,:,:,i) = tau_this;
    sigma2_all(:,:,:,i) = sigma2_this;

    CI_all(:,:,:,:,i) = CI_this;
    CIbc_all(:,:,:,:,i) = CIbc_this;

end

toc;


% ==================  Assemble results and output tables ==================

% load true T and sigma
file_name_true_T_sigma = ['Results/model_', num2str(model),'_true_T_sigma.mat'];
r_struct = load(file_name_true_T_sigma);
T_true = r_struct.True_T_table;
sigma2_true = r_struct.sigma2_est_table;

selected_pq = [1,2,4,5];
s_pq_num = length(selected_pq);
selected_n = 1:3;
s_n_num = length(selected_n);

cover_prob_all = zeros(s_pq_num, rho_num,s_n_num);
coverbc_prob_all = zeros(s_pq_num, rho_num, s_n_num);
RMSE_all = zeros(s_pq_num, rho_num, s_n_num);

for j0 = 1:s_pq_num

    j = selected_pq(j0);

    for k = 1:rho_num

        for r0 = 1: s_n_num
            r = selected_n(r0);

            T = T_true(j,k);

            sigma2 = sigma2_true(j,k);

            c_low = squeeze(CI_all(1,j,k,r,:));

            c_up = squeeze(CI_all(2,j,k,r,:));

            cover_indicator = c_low < T & c_up > T;

            cpb = mean(cover_indicator);

            cbc_low = squeeze(CIbc_all(1,j,k,r,:));

            cbc_up = squeeze(CIbc_all(2,j,k,r,:));

            cover_indicator_bc = cbc_low < T & cbc_up > T;

            cpb_bc = mean(cover_indicator_bc);

            cover_prob_all(j0,k,r0) = cpb;
            coverbc_prob_all(j0,k,r0) = cpb_bc;

            %-----
            sigma2_est_all = squeeze(sigma2_all(j,k,r,:));
            sser = (sigma2_est_all-sigma2)/sigma2;

            RMSE = sqrt(mean(sser.^2));

            RMSE_all(j0,k,r0) = RMSE;

        end

    end
end

Re_tensor = zeros(s_n_num, 3, rho_num, s_pq_num);

for i = 1:s_pq_num
    for j = 1:rho_num
        Re_tensor(:,1,j,i) = squeeze(cover_prob_all(i,j,:));
        Re_tensor(:,2,j,i) = squeeze(coverbc_prob_all(i,j,:));
        Re_tensor(:,3,j,i) = squeeze(RMSE_all(i,j,:));
    end
end

B = reshape(permute(Re_tensor,[1,3,2,4]), s_n_num * rho_num, 3 * s_pq_num);

save_table_name = ['Results/model_', num2str(model), '_rep_', num2str(rep_num),'_ECP_result_table.xlsx'];

C = compose('%.2f',B);
writecell(C, save_table_name)

save_all_name = ['Results/model_', num2str(model), '_rep_', num2str(rep_num),'_ECP_result.mat'];

save(save_all_name);

return;
