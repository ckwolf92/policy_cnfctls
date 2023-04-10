%% RUN POLICY COUNTERFACTUALS FOR INVESTMENT SHOCK
% Alisdair McKay & Christian Wolf
% this version: 03/24/2023

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

local = '/Users/christianwolf/Documents';
path = [local '/GitHub/policy_cnfctls/empirics'];
task = '/policy_counterfactuals';

addpath([path '/_auxiliary_functions'])
addpath([path '/_data'])
addpath([path task '/_aux'])
addpath([path '/var_inputs/_results'])

cd([path task]);

%% SETTINGS & IMPORTS

% for Figures 3-4 and C2-C3: use baseline monetary shock and cycle through
% counterfactual options

% for Figure C4: use alternative monetary shock and cycle through the
% output gap targeting and Taylor rule options

%----------------------------------------------------------------
% Import Investment Shocks
%----------------------------------------------------------------

load IS_inv

y_z  = -IS_inv.IRF_OLS(:,1);
pi_z = -IS_inv.IRF_OLS(:,2);
ib_z = -IS_inv.IRF_OLS(:,3);

T = length(y_z);

clear IS_inv

%----------------------------------------------------------------
% Select Monetary Shocks
%----------------------------------------------------------------

use_base = 1; % baseline monetary shocks (Romer-Romer, Gertler-Karadi)
use_alt  = 0; % alternative monetary shocks (Aruoba-Drechsel, MAR)

if use_base == 1
    load IS_mp_base
elseif use_alt == 1
    load IS_mp_alt
end

n_draws = size(IS_1.IRF,3);

%----------------------------------------------------------------
% Select Counterfactuals
%----------------------------------------------------------------

cnfctl_ngdp = 0; % nominal GDP targeting
cnfctl_0y   = 1; % output gap targeting
cnfctl_0ib  = 0; % nominal rate peg
cnfctl_tylr = 0; % Taylor rule

set_cnfctl_rule

%% CONSTRUCT COUNTERFACTUALS

disp('I am constructing the robust policy rule counterfactual.')

get_cnfctl

disp('Done!')

%% PLOT RESULTS

%----------------------------------------------------------------
% General Settings
%----------------------------------------------------------------

% plot horizon

IRF_hor_plot = 30;

% color settings

settings.colors.black  = [0 0 0];
settings.colors.dgrey  = [130/255 130/255 130/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];

% results scaling

scale = 5;

% go to results folder

cd([path task '/_results']);

%----------------------------------------------------------------
% Figure
%----------------------------------------------------------------

plotwidth = 0.27;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];

figure(1)

subplot(1,3,1)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(1:1:IRF_hor_plot,(scale*y_cnfctl_lb(1:IRF_hor_plot))',(scale*y_cnfctl_ub(1:IRF_hor_plot))',...
    settings.colors.lblue,settings.colors.lblue,0,1);
hold on
plot(1:1:IRF_hor_plot,scale*y_z(1:IRF_hor_plot),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*y_cnfctl_med(1:IRF_hor_plot),'Color',settings.colors.blue,'LineWidth',4)
hold on
if cnfctl_0y == 1
   plot(1:1:IRF_hor_plot,zeros(IRF_hor_plot,1),'--','Color',settings.colors.black,'LineWidth',4)
end
xlim([1 IRF_hor_plot])
ylim([-2 2])
yticks([-2 -1 0 1 2])
set(gcf,'color','w')
title('Output Gap','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
ylabel('\% Deviation','interpreter','latex','FontSize',20)
grid on
hold off

subplot(1,3,2)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(1:1:IRF_hor_plot,(scale*pi_cnfctl_lb(1:IRF_hor_plot))',(scale*pi_cnfctl_ub(1:IRF_hor_plot))',...
    settings.colors.lblue,settings.colors.lblue,0,1);
hold on
plot(1:1:IRF_hor_plot,scale*pi_z(1:IRF_hor_plot),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*pi_cnfctl_med(1:IRF_hor_plot),'Color',settings.colors.blue,'LineWidth',4)
hold on
if cnfctl_ngdp == 1
    plot(1:1:IRF_hor_plot,scale*pi_cnfctl_target_med(1:IRF_hor_plot),'--','Color',settings.colors.black,'LineWidth',4)
    hold on
    plot(1:1:IRF_hor_plot,scale*pi_cnfctl_target_lb(1:IRF_hor_plot),'--','Color',settings.colors.black,'LineWidth',2)
    hold on
    plot(1:1:IRF_hor_plot,scale*pi_cnfctl_target_ub(1:IRF_hor_plot),'--','Color',settings.colors.black,'LineWidth',2)
    hold on
end
xlim([1 IRF_hor_plot])
ylim([-2 2])
yticks([-2 -1 0 1 2])
set(gcf,'color','w')
title('Inflation','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
grid on
hold off

subplot(1,3,3)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
h = jbfill(1:1:IRF_hor_plot,(scale*ib_cnfctl_lb(1:IRF_hor_plot))',(scale*ib_cnfctl_ub(1:IRF_hor_plot))',...
    settings.colors.lblue,settings.colors.lblue,0,1);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on
plot(1:1:IRF_hor_plot,scale*ib_z(1:IRF_hor_plot),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*ib_cnfctl_med(1:IRF_hor_plot),'Color',settings.colors.blue,'LineWidth',4)
hold on
if cnfctl_tylr == 1
    plot(1:1:IRF_hor_plot,scale*ib_cnfctl_target_med(1:IRF_hor_plot),'--','Color',settings.colors.black,'LineWidth',4)
    hold on
    plot(1:1:IRF_hor_plot,scale*ib_cnfctl_target_lb(1:IRF_hor_plot),'--','Color',settings.colors.black,'LineWidth',2)
    hold on
    plot(1:1:IRF_hor_plot,scale*ib_cnfctl_target_ub(1:IRF_hor_plot),'--','Color',settings.colors.black,'LineWidth',2)
    hold on
elseif cnfctl_0y == 1
   plot(0,0,'--','Color',settings.colors.black,'LineWidth',4)
elseif cnfctl_0ib == 1
    plot(1:1:IRF_hor_plot,zeros(IRF_hor_plot,1),'--','Color',settings.colors.black,'LineWidth',4)
elseif cnfctl_ngdp == 1
   plot(0,0,'--','Color',settings.colors.black,'LineWidth',4)
end
xlim([1 IRF_hor_plot])
ylim([-3 1])
yticks([-3 -2 -1 0 1])
set(gcf,'color','w')
title('Interest Rate','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
if cnfctl_tylr == 1
    legend({'Baseline','Robust Counterfactual','Exact Taylor Rule'},'Location','Southeast','fontsize',18,'interpreter','latex')
elseif cnfctl_0y == 1
    legend({'Baseline','Robust Counterfactual','Exact Output Gap Targeting'},'Location','Southeast','fontsize',18,'interpreter','latex')
elseif cnfctl_0ib == 1
    legend({'Baseline','Robust Counterfactual','Exact Interest Rate Peg'},'Location','Southeast','fontsize',18,'interpreter','latex')
elseif cnfctl_ngdp == 1
    legend({'Baseline','Robust Counterfactual','Exact NGDP Targeting'},'Location','Southeast','fontsize',18,'interpreter','latex')
end
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if use_base == 1
    if cnfctl_0y == 1
        print('figure_3','-dpng');
    end
    if cnfctl_tylr == 1
        print('figure_4','-dpng');
    end
    if cnfctl_0ib == 1
        print('figure_c2','-dpng');
    end
    if cnfctl_ngdp == 1
        print('figure_c3','-dpng');
    end
elseif use_alt == 1
    if cnfctl_0y == 1
        print('figure_c4_1','-dpng');
    end
    if cnfctl_tylr == 1
        print('figure_c4_2','-dpng');
    end
end

cd([path task]);