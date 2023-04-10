clearvars -except path experiment

%----------------------------------------------------------------
% HANK Inputs
%----------------------------------------------------------------

load inputs_hank

%----------------------------------------------------------------
% Other Parameters
%----------------------------------------------------------------

% price stickiness

kappa = 0.1;

% government policy

phi_i_cnfctl  = 0.9;
phi_pi_cnfctl = 2;
phi_y_cnfctl  = 0.5;

phi_i_base  = 0;
phi_pi_base = 1.5;
phi_y_base  = 0;

%----------------------------------------------------------------
% Time Horizon
%----------------------------------------------------------------

T = size(C_ib,1);

%----------------------------------------------------------------
% Shock
%----------------------------------------------------------------

shock_s = zeros(T,1);

shock_s(1) = 1;
for t = 2:T
    shock_s(t) = 0.8 * shock_s(t-1);
end

cd([path experiment '/_results'])

save model_param -regexp ^(?!(path|experiment)$).

cd([path experiment]);