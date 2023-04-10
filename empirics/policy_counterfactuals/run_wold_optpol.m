%% RUN OPTIMAL POLICY FOR WOLD SHOCKS
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

% for Table C.1: use baseline monetary shock

%----------------------------------------------------------------
% Import Wold Shocks
%----------------------------------------------------------------

load IS_wold

Theta_base = IS_wold.Theta_OLS;
cov_base   = IS_wold.cov;
corr_base  = IS_wold.corr;

T   = size(Theta_base,3);
n_y = size(Theta_base,1);

clear IS_wold

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
% Set Policymaker Preferences
%----------------------------------------------------------------

set_polpref

%% CONSTRUCT COUNTERFACTUALS

disp('I am constructing the robust optimal policy counterfactual.')

get_optpol_wold

disp('Done!')

%% REPORT RESULTS

cd([path task '/_results']);

%----------------------------------------------------------------
% Baseline
%----------------------------------------------------------------

vars = {'y';'pi';'i'};
sd_base = sqrt(diag(cov_base));
corr_base = corr_base(:,1);

table_baseline = table(sd_base,corr_base,'RowNames',vars);
disp(table_baseline)

writetable(table_baseline)

%----------------------------------------------------------------
% Counterfactual
%----------------------------------------------------------------

vars = {'y';'pi';'i'};
sd_cnfctl = {[num2str(sqrt(cov_cnfctl_med(1,1))) ' (' num2str(sqrt(cov_cnfctl_lb(1,1))) ' ' num2str(sqrt(cov_cnfctl_ub(1,1))) ')'];...
    [num2str(sqrt(cov_cnfctl_med(2,2))) ' (' num2str(sqrt(cov_cnfctl_lb(2,2))) ' ' num2str(sqrt(cov_cnfctl_ub(2,2))) ')'];...
    [num2str(sqrt(cov_cnfctl_med(3,3))) ' (' num2str(sqrt(cov_cnfctl_lb(3,3))) ' ' num2str(sqrt(cov_cnfctl_ub(3,3))) ')']};
corr_cnfctl = {[num2str(corr_cnfctl_med(1,1)) ' (' num2str(corr_cnfctl_lb(1,1)) ' ' num2str(corr_cnfctl_ub(1,1)) ')'];...
    [num2str(corr_cnfctl_med(2,1)) ' (' num2str(corr_cnfctl_lb(2,1)) ' ' num2str(corr_cnfctl_ub(2,1)) ')'];...
    [num2str(corr_cnfctl_med(3,1)) ' (' num2str(corr_cnfctl_lb(3,1)) ' ' num2str(corr_cnfctl_ub(3,1)) ')']};

table_cnfctl = table(sd_cnfctl,corr_cnfctl,'RowNames',vars);
disp(table_cnfctl)

writetable(table_cnfctl)

cd([path task]);