function [] = plot_sigma_performance(model, true_sigma2_vec, sigma_est_all, tau_all, sigma_time_all, save_file_name)

selected_pq = [1,2,4,5];
spq_num = length(selected_pq);

n_all = [1000,2000,5000,10000];
n_num = length(n_all);

rRSME1_all = zeros(spq_num, n_num);
rRSME2_all = zeros(spq_num, n_num);

MT1_all = zeros(spq_num, n_num);
MT2_all = zeros(spq_num, n_num);

ECP1_all = zeros(spq_num, n_num);
ECP2_all = zeros(spq_num, n_num);

z_alpha_2 = norminv(0.975);


pqlabels = {'$(p,q)=(1,1)$','$(p,q)=(1,2)$','$(p,q)=(3,1)$','$(p,q)=(3,3)$','$(p,q)=(5,5)$'};


for pq_index = selected_pq

    true_sigma2 = true_sigma2_vec(pq_index);

    for n_index = 1:4

        this_n = n_all(n_index);


        ss = squeeze(sigma_est_all(:,pq_index,1,n_index,:));

        zz = squeeze(tau_all(:,pq_index,1,n_index,1));
        %zz_bc = squeeze(zeta_all(:,pq_index,1,n_index,2));


        ECP1_s = zz > -z_alpha_2 * ss(:,1)/sqrt(this_n) & zz < z_alpha_2 * ss(:,1)/sqrt(this_n);

        ECP2_s = zz > -z_alpha_2 * ss(:,2)/sqrt(this_n) & zz < z_alpha_2 * ss(:,2)/sqrt(this_n);

        ECP1 = mean(ECP1_s);

        ECP2 = mean(ECP2_s);

        ECP1_all(pq_index, n_index) = ECP1;

        ECP2_all(pq_index, n_index) = ECP2;


        ss2 = ss.^2;

        sser = (ss2 - true_sigma2)/true_sigma2;

        rRSME1 = sqrt(mean(sser(:,1).^2));

        rRSME2 = sqrt(mean(sser(:,2).^2));

        rRSME1_all(pq_index, n_index) = rRSME1;

        rRSME2_all(pq_index, n_index) = rRSME2;

        tt = squeeze(sigma_time_all(:,pq_index,1,n_index,:));

        MT1 = mean(tt(:,1));
        MT2 = mean(tt(:,2));

        MT1_all(pq_index, n_index) = MT1;
        MT2_all(pq_index, n_index) = MT2;



    end

end

lw = 1;
%ms = 5;

hAll = [];

figure('Units','inches','Position',[1,1,10,6]);

for i = 1:spq_num

    j = selected_pq(i);
    subplot(3,spq_num,i);
    h1 = plot(n_all, rRSME1_all(j,:),'b-*','LineWidth',lw,'MarkerSize',4);
    hold on;
    h2 = plot(n_all, rRSME2_all(j,:),'r-.o','LineWidth',lw,'MarkerSize',4);

    xlim([0,10200])
    xticks([1000,5000,10000])
    yticks(0:0.1:0.5);
    ylim([0,0.53]);
    xlabel('$n$','Interpreter','latex');
    ylabel('rRMSE');

    if j == 1
        legend({'$\hat{\sigma}_{0,\mathrm{F}}^2$','$\hat{\sigma}_{0,\mathrm{B}}^2$'}, 'Interpreter','latex','FontSize',10);
    end

    title( pqlabels{j},'Interpreter','latex');

    if isempty(hAll)
        hAll = [h1,h2];
    end

    %set(gca, 'FontSize',8)

    subplot(3,spq_num,i+spq_num);

    plot(n_all, ECP1_all(j,:),'b-*','LineWidth',lw,'MarkerSize',4);

    hold on;


    plot(n_all, ECP2_all(j,:),'r-.o','LineWidth',lw,'MarkerSize',4);

    %      if j < 4
    %      ylim([0.6,1])
    %
    %      yticks([0.6,0.95,1]);
    %      else
    %      ylim([0.3,1])
    %      yticks([0.3,0.6,1]);
    %
    %      end
    if model == 1
        ylim([0.25,1])
        yticks([0.3:0.2:0.9,1])
    else
        ylim([0,1])
        yticks([0,0.1:0.2:0.9,1])
    end


    xlim([0,10200])
    xticks([1000,5000,10000])

    xlabel('$n$','Interpreter','latex');
    ylabel('ECP');

    title( pqlabels{j},'Interpreter','latex');


    if j == 1
        legend({'$\hat{\sigma}_{0,\mathrm{F}}^2$','$\hat{\sigma}_{0,\mathrm{B}}^2$'},...
            'Interpreter','latex','FontSize',10,'Location','east');
    end




    subplot(3,spq_num,i+spq_num*2);

    %plot(n_all, MT1_all(j,:));
    plot(n_all, MT1_all(j,:),'b-*','LineWidth',lw,'MarkerSize',4);
    hold on;
    %plot(n_all, MT2_all(j,:));
    plot(n_all, MT2_all(j,:),'r-.o','LineWidth',lw,'MarkerSize',4);
    ylim([0,0.38]);
    yticks(0:0.1:0.3);

    xlim([0,10200]);
    xticks([1000,5000,10000])
    xlabel('$n$','Interpreter','latex');
    ylabel('CPU time (s)');

    if j == 1
        legend({'$\hat{\sigma}_{0,\mathrm{F}}^2$','$\hat{\sigma}_{0,\mathrm{B}}^2$'}, 'Interpreter','latex','FontSize',10);
    end

    title( pqlabels{j},'Interpreter','latex')



end

%result_name = ['Results/Model',num2str(model),'_sigma.eps'];

print(save_file_name,'-depsc');

% axlegend = axes('Position',[0 0 1 1],'Visible','off');
%
%  lgd =  legend(axlegend,hAll,{'asf','asdf'},'Interpreter','latex','FontSize',16, 'NumColumns',2,'Location','southoutside');
%
%  pos = lgd.Position;
%
%   pos(2) = pos(2) -0.01;
%
%   lgd.Position = pos;

