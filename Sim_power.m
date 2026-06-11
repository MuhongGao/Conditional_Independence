% Name: Sim_power.m
%
% Function: Conduct simulations realted to empirical powers in Section 5.2 (Figures 1-4)
%
% Input: Input the model number (1 or 2) and replication num, following the prompt.
%
% Note: This code performs parallel computing. Default Corenum = 64.
%       User may change the "Corenum" according to the machines capacity.
% 
% Elasped time: about 3 hours for rep_num = 1000 (Corenum = 64). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('Functions/');

model = input('Input the model number (1 or 2): ');
rep_num = input('Input the replication number (100, 500, or 1000): ');

%=====================  Settings ==============================
bootstrap_num = 200;

z_alpha = norminv(0.95);


p_all = [1,1,3,3,5];
q_all = [1,2,1,3,5];
pq_num = length(p_all);
rho_all = 0:0.1:0.8;
rho_num = length(rho_all);

n_vec = [1000,2000,5000,10000]; 
n_num = length(n_vec);
n_all = ones(pq_num, 1) * n_vec;

Result_all = zeros(rep_num, pq_num, rho_num,n_num,4);
tau_all = zeros(rep_num, pq_num, rho_num,n_num,2);

sigma_est_all = zeros(rep_num,pq_num,rho_num,n_num,2);
sigma_time_all = zeros(rep_num,pq_num, rho_num,n_num,2);


CoreNum = 64;

if isempty(gcp('nocreate'))
    parpool(CoreNum);
end

t0 = tic;

parfor (i = 1: rep_num, CoreNum)

    if mod(i,100) == 1
        fprintf([' rep = ', num2str(i), '\n']);
    end

    seed = i+100;
    rng(seed, 'twister');

    p_allc = p_all;
    q_allc = q_all;
    rho_allc = rho_all;
    n_allc = n_all;

    Result_this = zeros(pq_num,rho_num,n_num,4);
    tau_this = zeros(pq_num,rho_num,n_num,2);
    sigma_est_this = zeros(pq_num,rho_num,n_num,2);
    sigma_time_this = zeros(pq_num,rho_num,n_num,2);


    for pq_index = 1:pq_num
        
        for rho_index = 1:rho_num
            
            for n_index = 1:n_num

                p = p_allc(pq_index);
                q = q_allc(pq_index);
                rho = rho_allc(rho_index);
                n = n_allc(pq_index, n_index);

                m = ceil(n^0.5);

                lambda = n^(-0.85);

                if rho == 0
                    H = 0;
                else
                    H = 1;
                end


                [Z,Y,X] = data_gener(n,q,p, model, H, rho);

                %[tau, ~] = compute_zeta(Z, Y, X);

                [~, ~, ~, tau, ~] = compute_T(Z, Y, X, 0);

                Lhat1 = bias_correct([Z,X],Y, lambda);

                Lhat2 = bias_correct(Z,Y, lambda);

                bias = Lhat1- Lhat2; % do not multiply 6! tau is not xi, it is not normalized.

                tau_bc = tau - bias;

                % Limiting Variance estimation ===================
                %
                % method 1: direct estimation ------------
                t = tic;
                [sigma2_est_F] = est_sigma2_F(Z, Y, q, p);

                time_1 = toc(t);

                if sigma2_est_F >0
                    sigma_est_F = sqrt(sigma2_est_F);
                else
                    sigma_est_F = 0.01;
                end

                test_stat_F = sqrt(n) *tau / sigma_est_F;

                test_stat_F_bc = sqrt(n) *tau_bc / sigma_est_F;

                % method 2: m out of bootstrap -----------

                t = tic;
                [sigma2_est_B] = m_n_bootstrap(Z, Y, X, bootstrap_num, m);

                sigma_est_B = sqrt(sigma2_est_B);

                test_stat_B = sqrt(n) *tau / sigma_est_B;

                test_stat_B_bc =  sqrt(n) *tau_bc / sigma_est_B;

                time_2 = toc(t);

                %----- record decisions ------------
                Decision_this = zeros(1,4);

                Decision_this(1) = test_stat_F > z_alpha ;
                Decision_this(2) = test_stat_F_bc > z_alpha ;
                Decision_this(3) = test_stat_B > z_alpha ;
                Decision_this(4) = test_stat_B_bc > z_alpha ;

                %-------record tau ------------
                tau_this_0 = zeros(1,2);

                tau_this_0(1) = tau;
                tau_this_0(2) = tau_bc;



                Result_this(pq_index,rho_index,n_index, :) = reshape(Decision_this,1,1,1,4);
                tau_this(pq_index,rho_index,n_index, :) = reshape(tau_this_0,1,1,1,2);

                sigma_est_this(pq_index,rho_index,n_index, 1) = sigma_est_F;
                sigma_est_this(pq_index,rho_index,n_index, 2) = sigma_est_B;
                sigma_time_this(pq_index,rho_index,n_index, 1) = time_1;
                sigma_time_this(pq_index,rho_index,n_index, 2) = time_2;
            end
        end
    end

    Result_all(i,:,:,:,:) = reshape(Result_this, 1,pq_num,rho_num,n_num,4);
    tau_all(i,:,:,:,:) = reshape(tau_this, 1,pq_num,rho_num,n_num,2);
    sigma_est_all(i,:,:,:,:) = reshape(sigma_est_this, 1,pq_num,rho_num,n_num,2);
    sigma_time_all(i,:,:,:,:) = reshape(sigma_time_this, 1,pq_num,rho_num,n_num,2);

end



t_elapsed = toc(t0);

%==================== Output figures and tables =============================

Result_mean = mean(Result_all, 1);
R = squeeze(Result_mean);

TR = zeros(n_num*pq_num,4*rho_num);


for i = 1: pq_num
    for j = 1:rho_num
        R_this = squeeze(R(i,j,:,:));

        rows = n_num*(i-1) +1 :n_num*i;
        cols = 4*(j-1) +1 : 4*j;

        R_this = R_this(:,[1,3,2,4]);

        TR(rows, cols) = R_this;

    end
end

TG = compose('%.2f', TR);

% Plot power curves --------
save_file_name_1 = ['Results/Model',num2str(model),  '_rep_', num2str(rep_num), '_power.eps'];
plot_power_curves(R, save_file_name_1 );

% Plot rRMSE, ECP, CPU time, for comparing sigma estimators --------
file_name_true_T_sigma = ['Results/model_', num2str(model),'_true_T_sigma.mat'];
r_struct = load(file_name_true_T_sigma);
%T_true = r_struct.True_T_table;
sigma2_true = r_struct.sigma2_est_table;
true_sigma2_vec = sigma2_true(:,1);


save_file_name_2 = ['Results/Model',num2str(model), '_rep_', num2str(rep_num), '_sigma.eps'];
plot_sigma_performance(model, true_sigma2_vec, sigma_est_all, tau_all, sigma_time_all, save_file_name_2);

% Save all varibles
save_all_name = ['Results/model_', num2str(model), '_rep_', num2str(rep_num),'_power_result.mat'];

save(save_all_name);


return;


