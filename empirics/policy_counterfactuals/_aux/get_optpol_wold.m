%% COMPUTE COUNTERFACTUAL POLICY OUTCOMES

%----------------------------------------------------------------
% Counterfactual Wold IRFs
%----------------------------------------------------------------

Theta_cnfctl = NaN(n_y,n_y,T,n_draws);

for i_draw = 1:n_draws
    
% get IRFs

IS_pol.IRF = NaN(T,3,2);
IS_pol.IRF(:,:,1) = IS_1.IRF(:,:,i_draw);
IS_pol.IRF(:,:,2) = IS_2.IRF(:,:,i_draw);

Theta_y  = squeeze(IS_pol.IRF(:,1,:));
Theta_pi = squeeze(IS_pol.IRF(:,2,:));
Theta_ib = squeeze(IS_pol.IRF(:,3,:));

for i_y = 1:n_y

% get the corresponding base sequences

y_z  = squeeze(Theta_base(1,i_y,:));
pi_z = squeeze(Theta_base(2,i_y,:));
ib_z = squeeze(Theta_base(3,i_y,:));

% find best fit to counterfactual rule

[y_z_cnfctl,pi_z_cnfctl,ib_z_cnfctl] = optpol_fn(W_pi,W_y,...
    Theta_pi,Theta_y,Theta_ib,pi_z,y_z,ib_z);

Theta_cnfctl(1,i_y,:,i_draw) = y_z_cnfctl;
Theta_cnfctl(2,i_y,:,i_draw) = pi_z_cnfctl;
Theta_cnfctl(3,i_y,:,i_draw) = ib_z_cnfctl;

end

end

%----------------------------------------------------------------
% Counterfactual Second Moments
%----------------------------------------------------------------

% all posterior covariances and correlations

cov_cnfctl = zeros(n_y,n_y,n_draws);

for i_draw = 1:n_draws
    for i_hor = 1:T
        cov_cnfctl(:,:,i_draw) = cov_cnfctl(:,:,i_draw) + Theta_cnfctl(:,:,i_hor,i_draw) * Theta_cnfctl(:,:,i_hor,i_draw)';
    end
end

corr_cnfctl = zeros(n_y,n_y,n_draws);
for i_draw = 1:n_draws
    for i_y = 1:n_y
        for ii_y = 1:n_y
            corr_cnfctl(i_y,ii_y,i_draw) = cov_cnfctl(i_y,ii_y,i_draw)/sqrt(cov_cnfctl(i_y,i_y,i_draw) * cov_cnfctl(ii_y,ii_y,i_draw));
        end
    end
end

% percentiles

cov_cnfctl_lb  = zeros(n_y,n_y);
cov_cnfctl_med = zeros(n_y,n_y);
cov_cnfctl_ub  = zeros(n_y,n_y);

for i_y = 1:n_y
    for ii_y = 1:n_y
        cov_cnfctl_lb(i_y,ii_y)  = quantile(squeeze(cov_cnfctl(i_y,ii_y,:)),0.16);
        cov_cnfctl_med(i_y,ii_y) = quantile(squeeze(cov_cnfctl(i_y,ii_y,:)),0.5);
        cov_cnfctl_ub(i_y,ii_y)  = quantile(squeeze(cov_cnfctl(i_y,ii_y,:)),0.84);
    end
end

corr_cnfctl_lb  = zeros(n_y,n_y);
corr_cnfctl_med = zeros(n_y,n_y);
corr_cnfctl_ub  = zeros(n_y,n_y);

for i_y = 1:n_y
    for ii_y = 1:n_y
        corr_cnfctl_lb(i_y,ii_y)  = quantile(squeeze(corr_cnfctl(i_y,ii_y,:)),0.16);
        corr_cnfctl_med(i_y,ii_y) = quantile(squeeze(corr_cnfctl(i_y,ii_y,:)),0.5);
        corr_cnfctl_ub(i_y,ii_y)  = quantile(squeeze(corr_cnfctl(i_y,ii_y,:)),0.84);
    end
end