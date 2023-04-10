%% CONSTRUCT SUPPLY SHOCK COUNTERFACTUALS VIA MODEL SOLUTION VS. SUFFICIENT STATISTICS
% Alisdair McKay & Christian Wolf
% this version: 03/22/2023

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

local = '/Users/christianwolf/Documents';
path = [local '/GitHub/policy_cnfctls/theory'];
experiment = '/cnfctl_policy';

addpath([path '/_auxiliary_functions'])
addpath([path '/_hank_inputs'])
addpath([path experiment '/_aux'])
addpath([path experiment '/_results'])

cd([path experiment]);

%% MODEL PARAMETERS

disp('I am setting all model parameters.')

get_model_param

disp('Done!')

%% SOLVE MODEL UNDER COUNTERFACTUAL POLICY

disp('I am solving the model under the counterfactual policy rule.')

get_cnfctl_modelsolve

disp('Done!')

%% PREDICT COUNTERFACTUAL USING POLICY SHOCK CAUSAL EFFECTS UNDER BASELINE RULE

disp('I am using policy shock IRFs to predict the counterfactual.')

get_cnfctl_suffstats

disp('Done!')

%% PLOT RESULTS

%----------------------------------------------------------------
% Imports
%----------------------------------------------------------------

load cnfctl_modelsolve
load cnfctl_suffstats

%----------------------------------------------------------------
% Colors
%----------------------------------------------------------------

settings.colors.black = [0 0 0];
settings.colors.dgrey = [130/255 130/255 130/255];
settings.colors.blue  = [116/255 158/255 178/255];

settings.colors.var = zeros(40,3);
settings.colors.maxhor = 5;
settings.colors.maxweight = 0.8;
for i = 1:40
    if i <= settings.colors.maxhor
        weight = 0 + (i-1) * settings.colors.maxweight/(settings.colors.maxhor-1);
        settings.colors.var(i,:) = (1-weight) * settings.colors.black + weight * [1 1 1];
    else
        settings.colors.var(i,:) = (1-weight) * settings.colors.black + weight * [1 1 1];
    end
end
clear weight

%----------------------------------------------------------------
% Figure Settings
%----------------------------------------------------------------

% plot size

plotwidth = 0.25;
gapsize = 0.075;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * (gapsize + plotwidth)];

% IRF horizon

IRF_plot = 20;

% scale

scale = abs(1/y_s_base(1));

%----------------------------------------------------------------
% Plots
%----------------------------------------------------------------

cd([path experiment '/_results'])

for i_var_hor = 1:n_var_hor

figure(i_var_hor)

subplot(1,3,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:IRF_plot,scale*y_s_cnfctl(1:IRF_plot+1),'linewidth',5,'linestyle','-','color',settings.colors.blue)
hold on
plot(0:1:IRF_plot,scale*y_s_base(1:IRF_plot+1),'linewidth',5,'linestyle',':','color',settings.colors.dgrey)
hold on
plot(0:1:IRF_plot,scale*y_s_finite_all(1:IRF_plot+1,i_var_hor),'linewidth',5,'linestyle','--','color',settings.colors.black)
hold on
set(gcf,'color','w')
title('Output','interpreter','latex','fontsize',21)
xlabel('Horizon','interpreter','latex','FontSize',18)
ylabel('\%','interpreter','latex','FontSize',18,'Rotation',0)
ylim([-1 0.5])
grid on
hold off

subplot(1,3,2)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:IRF_plot,scale*pi_s_cnfctl(1:IRF_plot+1),'linewidth',5,'linestyle','-','color',settings.colors.blue)
hold on
plot(0:1:IRF_plot,scale*pi_s_base(1:IRF_plot+1),'linewidth',5,'linestyle',':','color',settings.colors.dgrey)
hold on
plot(0:1:IRF_plot,scale*pi_s_finite_all(1:IRF_plot+1,i_var_hor),'linewidth',5,'linestyle','--','color',settings.colors.black)
hold on
set(gcf,'color','w')
title('Inflation','interpreter','latex','fontsize',21)
xlabel('Horizon','interpreter','latex','FontSize',18)
legend({'True Counterfactual','Baseline','Predicted Counterfactual'},'Location','Northeast','fontsize',18,'interpreter','latex')
ylim([0 0.6])
grid on
hold off

subplot(1,3,3)
pos = get(gca, 'Position');
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
hold on
for i = 1:var_hor_all(i_var_hor)
    plot(0:1:IRF_plot,scale*shocks_s_finite_all(i,1:IRF_plot+1,i_var_hor),'linewidth',5,'linestyle','-','color',settings.colors.var(i,:))
    hold on
end
set(gcf,'color','w')
title('Policy Shocks','interpreter','latex','fontsize',21)
xlabel('Horizon','interpreter','latex','FontSize',18)
ylim([-0.6 0.1])
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.3*1.7*pos(3) 1.3*0.7*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print(strcat('figure_1_', num2str(i_var_hor)),'-dpng');

end

cd([path experiment]);