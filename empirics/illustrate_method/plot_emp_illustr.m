%% CONNECTION TO EMPIRICS: ILLUSTRATIVE FIGURE
% Alisdair McKay & Christian Wolf
% this version: 03/24/2023

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

local = '/Users/christianwolf/Documents';
path = [local '/GitHub/policy_cnfctls/empirics'];
task = '/illustrate_method';

cd([path task]);

%% INPUT

load emp_illustr_inputs.mat

%% FIGURES

%----------------------------------------------------------------
% Color Preparation
%----------------------------------------------------------------

settings.colors.black  = [0 0 0];
settings.colors.dgrey  = [130/255 130/255 130/255];
settings.colors.beige  = [196/255 174/255 120/255];
settings.colors.blue   = [116/255 158/255 178/255];

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

IRF_plot = 20;

%----------------------------------------------------------------
% Figure: Multiple Policy Shocks
%----------------------------------------------------------------

plotwidth = 0.18;
gapsize = 0.05;
gapsize_edges = (1-4*plotwidth-3*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * (gapsize + plotwidth), ...
    gapsize_edges + 3 * (gapsize + plotwidth)];

figure(1)

subplot(2,4,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:IRF_plot,pi_shock(1:IRF_plot+1),'linewidth',5,'linestyle','-','color',settings.colors.dgrey)
hold on
set(gcf,'color','w')
title('Oil Shock','interpreter','latex','fontsize',26)
ylabel('$\pi \; \;$','interpreter','latex','FontSize',22,'Rotation',0)
ylim([-0.5 1])
grid on
hold off

subplot(2,4,5)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:IRF_plot,i_shock(1:IRF_plot+1),'linewidth',5,'linestyle','-','color',settings.colors.dgrey)
hold on
set(gcf,'color','w')
xlabel('Horizon','interpreter','latex','FontSize',22)
ylabel('$i$','interpreter','latex','FontSize',22,'Rotation',0)
ylim([-0.2 0.4])
grid on
hold off

subplot(2,4,2)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:IRF_plot,pi_mp(1:IRF_plot+1),'linewidth',5,'linestyle','-','color',settings.colors.beige)
hold on
set(gcf,'color','w')
title('Policy Shock','interpreter','latex','fontsize',26)
ylim([-0.5 1])
grid on
hold off

subplot(2,4,6)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:IRF_plot,i_mp(1:IRF_plot+1),'linewidth',5,'linestyle','-','color',settings.colors.beige)
hold on
set(gcf,'color','w')
xlabel('Horizon','interpreter','latex','FontSize',22)
ylim([-0.2 0.4])
grid on
hold off

subplot(2,4,3)
pos = get(gca, 'Position');
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:IRF_plot,alpha*pi_mp_1(1:IRF_plot+1),'linewidth',5,'linestyle','-','color',settings.colors.beige)
hold on
plot(0:1:IRF_plot,pi_mp_2(1:IRF_plot+1),'linewidth',5,'linestyle','-.','color',settings.colors.beige)
hold on
set(gcf,'color','w')
title('Policy Shocks','interpreter','latex','fontsize',26)
ylim([-0.5 1])
grid on
hold off

subplot(2,4,7)
pos = get(gca, 'Position');
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:IRF_plot,alpha*i_mp_1(1:IRF_plot+1),'linewidth',5,'linestyle','-','color',settings.colors.beige)
hold on
plot(0:1:IRF_plot,i_mp_2(1:IRF_plot+1),'linewidth',5,'linestyle','-.','color',settings.colors.beige)
hold on
set(gcf,'color','w')
xlabel('Horizon','interpreter','latex','FontSize',22)
ylim([-0.2 0.4])
grid on
hold off

subplot(2,4,4)
pos = get(gca, 'Position');
pos(1) = left_pos(4);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:IRF_plot,pi_cnfctl(1:IRF_plot+1),'linewidth',5,'linestyle','-','color',settings.colors.blue)
hold on
set(gcf,'color','w')
title('Counterfactual','interpreter','latex','fontsize',26)
ylim([-0.5 1])
grid on
hold off

subplot(2,4,8)
pos = get(gca, 'Position');
pos(1) = left_pos(4);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:IRF_plot,i_cnfctl(1:IRF_plot+1),'linewidth',5,'linestyle','-','color',settings.colors.blue)
hold on
set(gcf,'color','w')
xlabel('Horizon','interpreter','latex','FontSize',22)
ylim([-0.2 0.4])
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.75*pos(3) 1.25*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_2','-dpng');