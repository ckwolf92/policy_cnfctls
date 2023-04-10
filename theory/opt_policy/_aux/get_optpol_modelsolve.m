%% CONSTRUCT SUPPLY SHOCK COUNTERFACTUALS VIA MODEL SOLUTION
% Alisdair McKay & Christian Wolf
% this version: 03/22/2023

%% HOUSEKEEPING

clearvars -except path experiment

%% PREPARATIONS

%----------------------------------------------------------------
% Global Variables
%----------------------------------------------------------------

% demand block

global C_ib C_pi C_y C_tau ...
    B_ib B_pi B_y B_tau ...
    beta r_b_SS y_tax Y_SS C_SS B_SS Trans_SS

% supply block

global kappa

% policy

global rho_b W Lambda

% shock path

global shock_s

% time horizon

global T

%----------------------------------------------------------------
% Imports
%----------------------------------------------------------------

load inputs_HANK
load model_param

%% CONSTRUCT PE JACOBIANS

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

step = 1;

%----------------------------------------------------------------
% NK Model Block
%----------------------------------------------------------------

% matrices used later for optimal policy FOCs

Pi_pi = zeros(T,T);
for t = 1:T
    Pi_pi(t,t) = 1;
    if t < T
        Pi_pi(t,t+1) = -1/(1+r_b_SS);
    end
end

Pi_y = kappa * eye(T);

% output

dpi_dy  = NaN(T,T);

for i_deriv = 1:T
    y_seq_deriv   = zeros(T,1);
    shock_s_deriv = zeros(T,1);
    y_seq_deriv(i_deriv,1) = y_seq_deriv(i_deriv,1) + step;
    pi_up = NKPC_fn(y_seq_deriv,shock_s_deriv);
    dpi_dy(:,i_deriv) = pi_up/step;
end

% supply shock

dpi_des = NaN(T,T);

for i_deriv = 1:T
    y_seq_deriv   = zeros(T,1);
    shock_s_deriv = zeros(T,1);
    shock_s_deriv(i_deriv,1) = shock_s_deriv(i_deriv,1) + step;
    pi_up = NKPC_fn(y_seq_deriv,shock_s_deriv);
    dpi_des(:,i_deriv) = pi_up/step;
end

%----------------------------------------------------------------
% Policy Rules
%----------------------------------------------------------------

dbg_dy   = NaN(T,T);
dtau_dy  = NaN(T,T);

for i_deriv = 1:T
    ib_seq_deriv   = zeros(T,1);
    pi_seq_deriv   = zeros(T,1);
    y_seq_deriv    = zeros(T,1);
    y_seq_deriv(i_deriv,1) = y_seq_deriv(i_deriv,1) + step;
    [taue_up,bg_up] = gov_bc_fn(ib_seq_deriv,pi_seq_deriv,y_seq_deriv);
    dbg_dy(:,i_deriv)  = bg_up/step;
    dtau_dy(:,i_deriv) = taue_up/step;
end

dbg_dpi  = NaN(T,T);
dtau_dpi = NaN(T,T);

for i_deriv = 1:T
    ib_seq_deriv   = zeros(T,1);
    pi_seq_deriv   = zeros(T,1);
    y_seq_deriv    = zeros(T,1);
    pi_seq_deriv(i_deriv,1) = pi_seq_deriv(i_deriv,1) + step;
    [taue_up,bg_up] = gov_bc_fn(ib_seq_deriv,pi_seq_deriv,y_seq_deriv);
    dbg_dpi(:,i_deriv)  = bg_up/step;
    dtau_dpi(:,i_deriv) = taue_up/step;
end

dbg_dib  = NaN(T,T);
dtau_dib = NaN(T,T);

for i_deriv = 1:T
    ib_seq_deriv   = zeros(T,1);
    pi_seq_deriv   = zeros(T,1);
    y_seq_deriv    = zeros(T,1);
    ib_seq_deriv(i_deriv,1) = ib_seq_deriv(i_deriv,1) + step;
    [taue_up,bg_up] = gov_bc_fn(ib_seq_deriv,pi_seq_deriv,y_seq_deriv);
    dbg_dib(:,i_deriv)  = bg_up/step;
    dtau_dib(:,i_deriv) = taue_up/step;
end

%----------------------------------------------------------------
% Consumers (for consistent notation)
%----------------------------------------------------------------

% nominal rate

dc_dib = C_ib;
db_dib = B_ib;

% inflation

dc_dpi = C_pi;
db_dpi = B_pi;

% income

dc_dy = C_y;
db_dy = B_y;

% transfer

dc_dtau = C_tau;
db_dtau = B_tau;

%% OPTIMAL POLICY RULE

%----------------------------------------------------------------
% FOC
%----------------------------------------------------------------

A = [lambda_pi * W, zeros(T,T), zeros(T,T), zeros(T,T), Pi_pi' * W, - C_pi' * W, - dtau_dpi' * W; ... % pi
    zeros(T,T), lambda_y * W, zeros(T,T), zeros(T,T), - Pi_y' * W, (eye(T) - C_y)' * W, - dtau_dy' * W; ... % y
    zeros(T,T), zeros(T,T), zeros(T,T), zeros(T,T), zeros(T,T), - C_tau' * W, W; ... % tau
    zeros(T,T), zeros(T,T), zeros(T,T), zeros(T,T), zeros(T,T), - C_ib' * W, -dtau_dib' * W]; % ib

%----------------------------------------------------------------
% Solve for Rule
%----------------------------------------------------------------

% use equations 2-4 to eliminate multipliers

A1_aux = -A(T+1:4*T,4*T+1:end);
A2_aux = A(T+1:4*T,1:4*T);
A3_aux = A(1:T,4*T+1:end) * A1_aux^(-1) * A2_aux;

% plug into first FOC

A_pi   = A(1:T,1:T) + A3_aux(:,1:T);
A_y    = A(1:T,T+1:2*T) + A3_aux(:,T+1:2*T);
A_tau  = A(1:T,2*T+1:3*T) + A3_aux(:,2*T+1:3*T);
A_ib   = A(1:T,3*T+1:4*T) + A3_aux(:,3*T+1:4*T);

% re-scale to get unit weights on inflation in policy criterion

scale = A_pi^(-1);

A_pi   = scale * A_pi;
A_y    = scale * A_y;
A_tau  = scale * A_tau;
A_ib   = scale * A_ib;

%% COMPUTE OPTIMAL POLICY IRFs

%----------------------------------------------------------------
% GE Updating Matrix
%----------------------------------------------------------------

M = [dbg_dy + dbg_dpi * dpi_dy - db_dy - db_dpi * dpi_dy - db_dtau * (dtau_dy + dtau_dpi * dpi_dy), dbg_dib - db_dib - db_dtau * dtau_dib; ...
    A_y + A_pi * dpi_dy + A_tau * (dtau_dy + dtau_dpi * dpi_dy), A_ib + A_tau * dtau_dib];

%----------------------------------------------------------------
% Get IRFs
%----------------------------------------------------------------

wedge_s = [(db_dpi * dpi_des + db_dtau * dtau_dpi * dpi_des - dbg_dpi * dpi_des) * shock_s; ...
            (-A_pi * dpi_des - A_tau * dtau_dpi * dpi_des) * shock_s];

sol_seq = M\wedge_s;

y_s    = sol_seq(1:T);
ib_s   = sol_seq(T+1:2*T);

pi_s    = dpi_dy * y_s + dpi_des * shock_s;
tau_s   = dtau_dpi * pi_s + dtau_dy * y_s + dtau_dib * ib_s;
bg_s    = dbg_dpi * pi_s + dbg_dy * y_s + dbg_dib * ib_s;
c_s     = dc_dy * y_s + dc_dpi * pi_s + dc_dtau * tau_s + dc_dib * ib_s;
b_s     = db_dy * y_s + db_dpi * pi_s + db_dtau * tau_s + db_dib * ib_s;
r_s     =  [0;ib_s(1:end-1)] - pi_s;

%% COLLECT SUPPLY SHOCK IRFs

y_s_optpol  = y_s;
pi_s_optpol = pi_s;
ib_s_optpol = ib_s;

A_pi_optpol = A_pi;
A_y_optpol  = A_y;
A_ib_optpol = A_ib;

cd([path experiment '/_results'])

save optpol_modelsolve y_s_optpol pi_s_optpol ib_s_optpol A_pi_optpol A_y_optpol A_ib_optpol
 
clearvars -except path experiment
cd([path experiment]);