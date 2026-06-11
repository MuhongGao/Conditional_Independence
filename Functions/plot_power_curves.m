function [] = plot_power_curves(R, save_file_name)

% lw1 = 1;
% lw = 1.4;

figure('Units','inches','Position',[1,1,7,7]);
%k=0;
hAll = [];

t = tiledlayout(5,4,'TileSpacing','compact','Padding','compact');

t.Position = [0.065,0.09,0.89,0.87];

axMat = gobjects(5,4);


xlabels = {'$n=1000$','$n=2000$','$n=5000$','$n=10000$'};
ylabels = {'$(p,q)=(1,1)$','$(p,q)=(1,2)$','$(p,q)=(3,1)$','$(p,q)=(3,3)$','$(p,q)=(5,5)$'};

method_names = {'$\mathsf{T}^{\mathrm{F}}$ \quad','$\mathsf{T}^{\mathrm{B}}$ \quad',...
    '$\mathsf{T}^{\mathrm{F,bc}}$ \quad','$\mathsf{T}^{\mathrm{B,bc}}$'};


for j = 1:5
    for i = 1: 4

        ax = nexttile;
        axMat(j,i) = ax;


        %k = k+1;
        %subplot(5,4,k);
        R_t = squeeze(R(j,:,i,:));
        R_t = R_t(:,[1,3,2,4]);

        h1=plot(0:0.1:0.8, R_t(:,1),'g-.','LineWidth', 1.9);
        hold on;
        h2=plot(0:0.1:0.8, R_t(:,2),'r--','LineWidth', 1.2);
        hold on;
        h3=plot(0:0.1:0.8, R_t(:,3),'c-', 'LineWidth', 1.8);
        hold on;
        %plot(0:0.1:0.8, R_t(:,4),':','Color',[1,0.5,0],'LineWidth', lw);
        h4=plot(0:0.1:0.8, R_t(:,4),'k:','LineWidth', 1);

        if isempty(hAll)
            hAll = [h1,h2,h3,h4];
        end

        xlim([0,0.8])
        xticks(0:0.2:0.8)
        yticks(0:0.2:1)
        xlabel('$\rho$', 'Interpreter','latex')
        ylabel('power')
        set(gca, 'FontSize',6)


    end

    axleft = axMat(j,1);
    %t.Children((j-1)*4+1);
    xl=axleft.XLim;
    yl = axleft.YLim;

    text(xl(1)-4.17*(xl(2)-xl(1)),mean(yl),ylabels{j},...
        'Interpreter','latex','Rotation',90,...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10)

    %yl = ylabel(t,'sdf');
    %yl.Layout.Tile = 'west';
    %yl.FontSize = 12;
end


axlegend = axes('Position',[0 0 1 1],'Visible','off');

% text(0.5,0.5,'sdffsa')
%
% for c = 1:5
%     pos = axMat(c,1).Position;
%     text(pos(1)-0.05,pos(2)+pos(4)/2-0.03,'$(p,q)=(1,2)$',...
%         'Interpreter','latex','Rotation',90,...
%         'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);
% end

for c = 1:4
    pos = axMat(1,c).Position;
    text(pos(1)+pos(3)/2,pos(2) + pos(4) +0.02,xlabels{c},...
        'Interpreter','latex',...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);
    %
    % text(pos(1),pos(2),'sdfsf');
end



lgd =  legend(axlegend,hAll,method_names,'Interpreter','latex','FontSize',11, 'NumColumns',4,'Location','southoutside');

pos = lgd.Position;

pos(2) = pos(2) +0.01;

lgd.Position = pos;

%set(gcf,'Position',[100,100,700,700])

%result_name = ['Results/Model',num2str(model),'_power.eps'];

print(save_file_name,'-depsc');

end


