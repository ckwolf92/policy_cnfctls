%% COMPUTE COUNTERFACTUAL POLICY OUTCOMES

y_optpol_draws    = NaN(T,n_draws);
pi_optpol_draws   = NaN(T,n_draws);
ib_optpol_draws   = NaN(T,n_draws);

for i_draw = 1:n_draws
    
% get IRFs

IS_pol.IRF = NaN(T,3,2);
IS_pol.IRF(:,:,1) = IS_1.IRF(:,:,i_draw);
IS_pol.IRF(:,:,2) = IS_2.IRF(:,:,i_draw);

% find optimal outcome

Theta_y  = squeeze(IS_pol.IRF(:,1,:));
Theta_pi = squeeze(IS_pol.IRF(:,2,:));
Theta_ib = squeeze(IS_pol.IRF(:,3,:));

[y_z_var,pi_z_var,ib_z_var] = optpol_fn(W_pi,W_y,...
    Theta_pi,Theta_y,Theta_ib,pi_z,y_z,ib_z);

y_optpol_draws(:,i_draw)    = y_z_var;
pi_optpol_draws(:,i_draw)   = pi_z_var;
ib_optpol_draws(:,i_draw)   = ib_z_var;

end

% collect IRF results

y_optpol_lb     = NaN(T,1);
y_optpol_ub     = NaN(T,1);
y_optpol_med    = NaN(T,1);
pi_optpol_lb    = NaN(T,1);
pi_optpol_ub    = NaN(T,1);
pi_optpol_med   = NaN(T,1);
ib_optpol_lb    = NaN(T,1);
ib_optpol_ub    = NaN(T,1);
ib_optpol_med   = NaN(T,1);

for t = 1:T
    y_optpol_med(t)    = quantile(y_optpol_draws(t,:),0.5);
    y_optpol_lb(t)     = quantile(y_optpol_draws(t,:),0.16);
    y_optpol_ub(t)     = quantile(y_optpol_draws(t,:),0.84);
    pi_optpol_med(t)   = quantile(pi_optpol_draws(t,:),0.5);
    pi_optpol_lb(t)    = quantile(pi_optpol_draws(t,:),0.16);
    pi_optpol_ub(t)    = quantile(pi_optpol_draws(t,:),0.84);
    ib_optpol_med(t)   = quantile(ib_optpol_draws(t,:),0.5);
    ib_optpol_lb(t)    = quantile(ib_optpol_draws(t,:),0.16);
    ib_optpol_ub(t)    = quantile(ib_optpol_draws(t,:),0.84);
end