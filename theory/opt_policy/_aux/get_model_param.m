clearvars -except path experiment

%----------------------------------------------------------------
% HANK Inputs
%----------------------------------------------------------------

load inputs_hank

%----------------------------------------------------------------
% Time Horizon
%----------------------------------------------------------------

T = size(C_ib,1);

%----------------------------------------------------------------
% Other Parameters
%----------------------------------------------------------------

% price stickiness

kappa = 0.1;

% government policy

lambda_pi = 1;
lambda_y  = 1;

Lambda = diag([lambda_pi,lambda_y]);

W = zeros(T,T);
for t = 1:T
    W(t,t) = (1/(1+r_b_SS))^(t-1);
end

phi_i_base  = 0;
phi_pi_base = 1.5;
phi_y_base  = 0;

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