%% GET HANK INPUTS FOR POLICY COUNTERFACTUAL ILLUSTRATION
% Alisdair McKay & Christian Wolf
% this version: 03/22/2023

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

local = '/Users/christianwolf/Documents';
path = [local '/GitHub/policy_cnfctls/theory/_hank_inputs'];

addpath(genpath([path '/_aux']))
addpath([path '/_income_process'])
addpath([path '/_steady_state'])
addpath([path '/_jacobians'])

cd(path);

%% COMPUTE STEADY STATE

disp('I am solving for the steady state of the model.')

get_ss

disp('Done!')

%% COMPUTE JACOBIANS

disp('I am computing the PE Jacobians.')

get_jacobians

disp('Done!')