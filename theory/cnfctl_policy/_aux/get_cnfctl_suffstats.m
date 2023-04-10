%% CONSTRUCT SUPPLY SHOCK COUNTERFACTUALS VIA SUFFICIENT STATISTICS
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

global phi_i phi_pi phi_y

% shock path

global shock_s

% time horizon

global T

%----------------------------------------------------------------
% Imports
%----------------------------------------------------------------

load inputs_HANK
load model_param
load cnfctl_modelsolve

%----------------------------------------------------------------
% Set Policy
%----------------------------------------------------------------

phi_i  = phi_i_base;
phi_pi = phi_pi_base;
phi_y  = phi_y_base;

%% CONSTRUCT PE JACOBIANS

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

step = 1;

%----------------------------------------------------------------
% NK Model Block
%----------------------------------------------------------------

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

% monetary policy

dib_dpi = NaN(T,T);

for i_deriv = 1:T
    pi_seq_deriv = zeros(T,1);
    y_seq_deriv  = zeros(T,1);
    m_seq_deriv  = zeros(T,1);
    pi_seq_deriv(i_deriv,1) = pi_seq_deriv(i_deriv,1) + step;
    ib_up = ib_rule_fn(pi_seq_deriv,y_seq_deriv,m_seq_deriv);
    dib_dpi(:,i_deriv) = ib_up/step;
end

dib_dy = NaN(T,T);

for i_deriv = 1:T
    pi_seq_deriv = zeros(T,1);
    y_seq_deriv  = zeros(T,1);
    m_seq_deriv  = zeros(T,1);
    y_seq_deriv(i_deriv,1) = y_seq_deriv(i_deriv,1) + step;
    ib_up = ib_rule_fn(pi_seq_deriv,y_seq_deriv,m_seq_deriv);
    dib_dy(:,i_deriv) = ib_up/step;
end

dib_dm = NaN(T,T);

for i_deriv = 1:T
    pi_seq_deriv = zeros(T,1);
    y_seq_deriv  = zeros(T,1);
    m_seq_deriv  = zeros(T,1);
    m_seq_deriv(i_deriv,1) = m_seq_deriv(i_deriv,1) + step;
    ib_up = ib_rule_fn(pi_seq_deriv,y_seq_deriv,m_seq_deriv);
    dib_dm(:,i_deriv) = ib_up/step;
end

% fiscal policy

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

%% GE IRF COMPUTATION

%----------------------------------------------------------------
% GE Multiplier
%----------------------------------------------------------------

M = (dbg_dy - db_dy + (dbg_dpi - db_dpi) * dpi_dy + (dbg_dib - db_dib) * (dib_dpi * dpi_dy + dib_dy) - db_dtau * (dtau_dy + dtau_dpi * dpi_dy ...
    + dtau_dib * (dib_dpi * dpi_dy + dib_dy)));

%----------------------------------------------------------------
% Supply Shock IRFs
%----------------------------------------------------------------

wedge_s = ((db_dpi-dbg_dpi) * dpi_des + (db_dib-dbg_dib) * dib_dpi * dpi_des ...
    + db_dtau * dtau_dpi * dpi_des + db_dtau * dtau_dib * dib_dpi * dpi_des);

Y_s = M^(-1) * wedge_s;

Pi_s   = dpi_dy * Y_s + dpi_des;
Ib_s   = dib_dpi * Pi_s + dib_dy * Y_s;
Tau_s  = dtau_dy * Y_s + dtau_dpi * Pi_s + dtau_dib * Ib_s;
D_s    = zeros(T,T);

X_s = [Pi_s; Y_s];

%----------------------------------------------------------------
% Interest Rate Policy IRFs
%----------------------------------------------------------------

wedge_ib = (db_dib - dbg_dib + db_dtau * dtau_dib) * dib_dm;

Y_ib = M^(-1) * wedge_ib;

Pi_ib   = dpi_dy * Y_ib;
Ib_ib   = dib_dm + dib_dpi * Pi_ib + dib_dy * Y_ib;
Tau_ib  = dtau_dy * Y_ib + dtau_dpi * Pi_ib + dtau_dib * Ib_ib;
D_ib    = zeros(T,T);

X_ib = [Pi_ib; Y_ib];

%% COLLECT SUPPLY SHOCK IRFs

y_s_base  = Y_s * shock_s;
pi_s_base = Pi_s * shock_s;
ib_s_base = Ib_s * shock_s;

%% PREDICT COUNTERFACTUAL USING POLICY SHOCKS

%----------------------------------------------------------------
% Full Sufficient Statistics Result
%----------------------------------------------------------------

aux_mat = [eye(T), zeros(T,T), zeros(T,T), -Y_ib; ...
    zeros(T,T), eye(T), zeros(T,T), -Pi_ib; ...
    zeros(T,T), zeros(T,T), eye(T), -Ib_ib; ...
    A_y_cnfctl, A_pi_cnfctl, A_ib_cnfctl, zeros(T,T)];

base_irf_mat   = [Y_s; Pi_s; Ib_s; zeros(T,T)];
target_irf_mat = aux_mat^(-1) * base_irf_mat;

cnfctl_irf_pred = target_irf_mat * shock_s;

y_s_cnfctl_pred  = cnfctl_irf_pred(1:T);
pi_s_cnfctl_pred = cnfctl_irf_pred(T+1:2*T);
ib_s_cnfctl_pred = cnfctl_irf_pred(2*T+1:3*T);
shock_cnfctl     = cnfctl_irf_pred(3*T+1:4*T);

%----------------------------------------------------------------
% Finite-Shock Generalization of Sims-Zha
%----------------------------------------------------------------

var_hor_all = [1 2 4 8 20 40];
n_var_hor   = length(var_hor_all);

y_s_finite_all  = NaN(T,n_var_hor);
pi_s_finite_all = NaN(T,n_var_hor);
ib_s_finite_all = NaN(T,n_var_hor);

shocks_s_finite_all = NaN(T,T,n_var_hor);

for i_var_hor = 1:n_var_hor

    var_hor = var_hor_all(i_var_hor);

    [y_s_var,pi_s_var,ib_s_var,shocks_var] ...
            = cnfctl_finite(var_hor,A_ib_cnfctl,A_pi_cnfctl,A_y_cnfctl,...
            Pi_ib,Y_ib,Ib_ib,Pi_s,Y_s,Ib_s,shock_s);

    y_s_finite_all(:,i_var_hor)  = y_s_var;
    pi_s_finite_all(:,i_var_hor) = pi_s_var;
    ib_s_finite_all(:,i_var_hor) = ib_s_var;
    shocks_s_finite_all(:,:,i_var_hor) = shocks_var;

end

shocks_s_finite_all = shocks_s_finite_all(1:40,1:40,:);

%% SAVE RESULTS

cd([path experiment '/_results'])

save cnfctl_suffstats y_s_base pi_s_base ib_s_base y_s_cnfctl_pred pi_s_cnfctl_pred ib_s_cnfctl_pred shock_cnfctl ...
    y_s_finite_all pi_s_finite_all ib_s_finite_all shocks_s_finite_all var_hor_all n_var_hor

clearvars -except path experiment
cd([path experiment]);