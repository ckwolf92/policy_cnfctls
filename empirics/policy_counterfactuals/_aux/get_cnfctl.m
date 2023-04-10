%% COMPUTE COUNTERFACTUAL POLICY OUTCOMES

y_cnfctl_draws    = NaN(T,n_draws);
pi_cnfctl_draws   = NaN(T,n_draws);
ib_cnfctl_draws   = NaN(T,n_draws);
rule_error_draws  = NaN(T,n_draws);
nu_cnfctl_draws   = NaN(2,n_draws);

for i_draw = 1:n_draws
    
% get IRFs

IS_pol.IRF = NaN(T,3,2);
IS_pol.IRF(:,:,1) = IS_1.IRF(:,:,i_draw);
IS_pol.IRF(:,:,2) = IS_2.IRF(:,:,i_draw);

% find best fit to counterfactual rule

Theta_y  = squeeze(IS_pol.IRF(:,1,:));
Theta_pi = squeeze(IS_pol.IRF(:,2,:));
Theta_ib = squeeze(IS_pol.IRF(:,3,:));

[y_z_var,pi_z_var,ib_z_var,nu_z_var,rule_error_var] = cnfctl_fn(A_y,A_pi,A_ib,...
    Theta_pi,Theta_y,Theta_ib,pi_z,y_z,ib_z);

y_cnfctl_draws(:,i_draw)    = y_z_var;
pi_cnfctl_draws(:,i_draw)   = pi_z_var;
ib_cnfctl_draws(:,i_draw)   = ib_z_var;
rule_error_draws(:,i_draw)  = rule_error_var;
nu_cnfctl_draws(:,i_draw)   = nu_z_var;

end

% collect IRF results

y_cnfctl_lb     = NaN(T,1);
y_cnfctl_ub     = NaN(T,1);
y_cnfctl_med    = NaN(T,1);
pi_cnfctl_lb    = NaN(T,1);
pi_cnfctl_ub    = NaN(T,1);
pi_cnfctl_med   = NaN(T,1);
ib_cnfctl_lb    = NaN(T,1);
ib_cnfctl_ub    = NaN(T,1);
ib_cnfctl_med   = NaN(T,1);
rule_error_lb   = NaN(T,1);
rule_error_ub   = NaN(T,1);
rule_error_med  = NaN(T,1);

for t = 1:T
    y_cnfctl_med(t)    = quantile(y_cnfctl_draws(t,:),0.5);
    y_cnfctl_lb(t)     = quantile(y_cnfctl_draws(t,:),0.16);
    y_cnfctl_ub(t)     = quantile(y_cnfctl_draws(t,:),0.84);
    pi_cnfctl_med(t)   = quantile(pi_cnfctl_draws(t,:),0.5);
    pi_cnfctl_lb(t)    = quantile(pi_cnfctl_draws(t,:),0.16);
    pi_cnfctl_ub(t)    = quantile(pi_cnfctl_draws(t,:),0.84);
    ib_cnfctl_med(t)   = quantile(ib_cnfctl_draws(t,:),0.5);
    ib_cnfctl_lb(t)    = quantile(ib_cnfctl_draws(t,:),0.16);
    ib_cnfctl_ub(t)    = quantile(ib_cnfctl_draws(t,:),0.84);
    rule_error_med(t)  = quantile(rule_error_draws(t,:),0.5);
    rule_error_lb(t)   = quantile(rule_error_draws(t,:),0.16);
    rule_error_ub(t)   = quantile(rule_error_draws(t,:),0.84);
end

% construct rule targets (if applicable)

if cnfctl_tylr == 1
    ib_cnfctl_target = NaN(T,n_draws);
    for i_draw = 1:n_draws
        ib_cnfctl_target(:,i_draw) = -A_ib^(-1) * (A_y * y_cnfctl_draws(:,i_draw) + A_pi * pi_cnfctl_draws(:,i_draw));
    end
    ib_cnfctl_target_med = NaN(T,1);
    ib_cnfctl_target_lb  = NaN(T,1);
    ib_cnfctl_target_ub  = NaN(T,1);
    for t = 1:T
        ib_cnfctl_target_med(t) = quantile(ib_cnfctl_target(t,:),0.5);
        ib_cnfctl_target_lb(t)  = quantile(ib_cnfctl_target(t,:),0.16);
        ib_cnfctl_target_ub(t)  = quantile(ib_cnfctl_target(t,:),0.84);
    end
elseif cnfctl_ngdp == 1
    pi_cnfctl_target = NaN(T,n_draws);
    for i_draw = 1:n_draws
        pi_cnfctl_target(:,i_draw) = -(y_cnfctl_draws(:,i_draw) - [0;y_cnfctl_draws(1:end-1,i_draw)]);
    end
    pi_cnfctl_target_med = NaN(T,1);
    pi_cnfctl_target_lb  = NaN(T,1);
    pi_cnfctl_target_ub  = NaN(T,1);
    for t = 1:T
        pi_cnfctl_target_med(t) = quantile(pi_cnfctl_target(t,:),0.5);
        pi_cnfctl_target_lb(t)  = quantile(pi_cnfctl_target(t,:),0.16);
        pi_cnfctl_target_ub(t)  = quantile(pi_cnfctl_target(t,:),0.84);
    end
end