%% MONETARY POLICY SHOCK IRFs: PLOT RESULTS
% Alisdair McKay & Christian Wolf
% this version: 03/24/2023

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

local = '/Users/christianwolf/Documents';
path = [local '/GitHub/policy_cnfctls/empirics'];
task = '/var_inputs';

addpath([path '/_auxiliary_functions'])
addpath([path '/_data'])
addpath([path task '/_results'])

cd([path task]);

%% IMPORT RESULTS

load IS_mp_base

IS_gk = IS_1;
IS_rr = IS_2;

clear IS_1 IS_2

%% PLOT RESULTS

%----------------------------------------------------------------
% General Settings
%----------------------------------------------------------------

% plot horizon

IRF_hor_plot = 28;

% color settings

settings.colors.black    = [0 0 0];
settings.colors.grey     = [230/255 230/255 230/255];
settings.colors.navyblue = [0/255 0/255 50/255];

% variable names

series_names = {'Output Gap','Inflation','Interest Rate'};

% plot size

plotwidth = 0.25;
gapsize = 0.075;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * (gapsize + plotwidth)];

%----------------------------------------------------------------
% Romer-Romer Plot
%----------------------------------------------------------------

cd([path task '/_results'])

figure(1)

subplot(1,3,1)
pos = get(gca, 'Position');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(1:1:IRF_hor_plot,(squeeze(IS_rr.IRF_lb(1:IRF_hor_plot,1)))',...
    (squeeze(IS_rr.IRF_ub(1:IRF_hor_plot,1)))',settings.colors.grey,settings.colors.grey,0,1);
plot(1:1:IRF_hor_plot,IS_rr.IRF_med(1:IRF_hor_plot,1),'Color',settings.colors.navyblue,'LineWidth',6)
hold on
hold on
xlim([1 IRF_hor_plot])
ylim([-0.35 0.15])
yticks([-0.3 -0.2 -0.1 0 0.1])
set(gcf,'color','w')
title('Output Gap','interpreter','latex','fontsize',26)
xlabel('Horizon','interpreter','latex','FontSize',22)
ylabel('\% Deviation','interpreter','latex','FontSize',22)
grid on
hold off

subplot(1,3,2)
pos = get(gca, 'Position');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(1:1:IRF_hor_plot,(squeeze(IS_rr.IRF_lb(1:IRF_hor_plot,2)))',...
    (squeeze(IS_rr.IRF_ub(1:IRF_hor_plot,2)))',settings.colors.grey,settings.colors.grey,0,1);
plot(1:1:IRF_hor_plot,IS_rr.IRF_med(1:IRF_hor_plot,2),'Color',settings.colors.navyblue,'LineWidth',6)
hold on
hold on
xlim([1 IRF_hor_plot])
ylim([-0.35 0.1])
yticks([-0.3 -0.2 -0.1 0 0.1])
set(gcf,'color','w')
title('Inflation','interpreter','latex','fontsize',26)
xlabel('Horizon','interpreter','latex','FontSize',22)
grid on
hold off

subplot(1,3,3)
pos = get(gca, 'Position');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(1:1:IRF_hor_plot,(squeeze(IS_rr.IRF_lb(1:IRF_hor_plot,3)))',...
    (squeeze(IS_rr.IRF_ub(1:IRF_hor_plot,3)))',settings.colors.grey,settings.colors.grey,0,1);
plot(1:1:IRF_hor_plot,IS_rr.IRF_med(1:IRF_hor_plot,3),'Color',settings.colors.navyblue,'LineWidth',6)
hold on
hold on
xlim([1 IRF_hor_plot])
ylim([-0.4 0.8])
yticks([-0.4 -0.2 0 0.2 0.4 0.6 0.8])
set(gcf,'color','w')
title('Interest Rate','interpreter','latex','fontsize',26)
xlabel('Horizon','interpreter','latex','FontSize',22)
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.72*pos(3) 1.12*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_c1_1','-dpng');

%----------------------------------------------------------------
% Gertler-Karadi Plot
%----------------------------------------------------------------

figure(2)

subplot(1,3,1)
pos = get(gca, 'Position');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(1:1:IRF_hor_plot,(squeeze(IS_gk.IRF_lb(1:IRF_hor_plot,1)))',...
    (squeeze(IS_gk.IRF_ub(1:IRF_hor_plot,1)))',settings.colors.grey,settings.colors.grey,0,1);
plot(1:1:IRF_hor_plot,IS_gk.IRF_med(1:IRF_hor_plot,1),'Color',settings.colors.navyblue,'LineWidth',6)
hold on
hold on
xlim([1 IRF_hor_plot])
ylim([-0.2 0.2])
yticks([-0.2 -0.1 0 0.1 0.2])
set(gcf,'color','w')
title('Output Gap','interpreter','latex','fontsize',26)
xlabel('Horizon','interpreter','latex','FontSize',22)
ylabel('\% Deviation','interpreter','latex','FontSize',22)
grid on
hold off

subplot(1,3,2)
pos = get(gca, 'Position');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(1:1:IRF_hor_plot,(squeeze(IS_gk.IRF_lb(1:IRF_hor_plot,2)))',...
    (squeeze(IS_gk.IRF_ub(1:IRF_hor_plot,2)))',settings.colors.grey,settings.colors.grey,0,1);
plot(1:1:IRF_hor_plot,IS_gk.IRF_med(1:IRF_hor_plot,2),'Color',settings.colors.navyblue,'LineWidth',6)
hold on
hold on
xlim([1 IRF_hor_plot])
ylim([-0.25 0.2])
yticks([-0.2 -0.1 0 0.1 0.2])
set(gcf,'color','w')
title('Inflation','interpreter','latex','fontsize',26)
xlabel('Horizon','interpreter','latex','FontSize',22)
grid on
hold off

subplot(1,3,3)
pos = get(gca, 'Position');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(1:1:IRF_hor_plot,(squeeze(IS_gk.IRF_lb(1:IRF_hor_plot,3)))',...
    (squeeze(IS_gk.IRF_ub(1:IRF_hor_plot,3)))',settings.colors.grey,settings.colors.grey,0,1);
plot(1:1:IRF_hor_plot,IS_gk.IRF_med(1:IRF_hor_plot,3),'Color',settings.colors.navyblue,'LineWidth',6)
hold on
hold on
xlim([1 IRF_hor_plot])
ylim([-0.2 0.6])
yticks([-0.2 0 0.2 0.4 0.6])
set(gcf,'color','w')
title('Interest Rate','interpreter','latex','fontsize',26)
xlabel('Horizon','interpreter','latex','FontSize',22)
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.72*pos(3) 1.12*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_c1_2','-dpng');

cd([path task]);