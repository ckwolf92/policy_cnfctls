%% RUN VAR INPUTS FOR POLICY COUNTERFACTUAL ANALYSIS
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
addpath([path task '/_aux'])

cd([path task]);

%% BASELINE MONETARY POLICY VAR

disp('I am running the baseline monetary policy VAR.')

run_mp_base

disp('Done!')

%% ALTERNATIVE MONETARY POLICY VAR

disp('I am running the alternative monetary policy VAR.')

run_mp_alt

disp('Done!')

%% INVESTMENT SHOCK VAR

disp('I am running the investment shock VAR.')

run_inv

disp('Done!')

%% WOLD VAR

disp('I am running a VAR for Wold shock IRFs.')

run_wold

disp('Done!')