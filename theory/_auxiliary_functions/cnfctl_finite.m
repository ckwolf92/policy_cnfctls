function [y_s_pred,pi_s_pred,ib_s_pred,nu_all,error_s_all,ib_s_surprise] = cnfctl_finite(h,A_ib,A_pi,A_y,...
    Pi_ib_true,Y_ib_true,Ib_ib_true,Pi_s,Y_s,Ib_s,shock_s);

% settings and placeholders

T = size(A_ib,1);

nu_all = zeros(T,T);
y_s_all = NaN(T,T);
pi_s_all = NaN(T,T);
ib_s_all = NaN(T,T);
error_s_all = NaN(T,T);

% construct impact shock sequence and IRFs

nu_all(1:h,1) = -(A_pi(1:h,:) * Pi_ib_true(:,1:h) + A_y(1:h,:) * Y_ib_true(:,1:h) ...
    + A_ib(1:h,:) * Ib_ib_true(:,1:h))^(-1) * (A_pi(1:h,:) * (Pi_s * shock_s) ...
    + A_y(1:h,:) * (Y_s * shock_s) + A_ib(1:h,:) * (Ib_s * shock_s));

y_s_all(:,1)  = Y_s * shock_s + Y_ib_true * nu_all(:,1);
pi_s_all(:,1) = Pi_s * shock_s + Pi_ib_true * nu_all(:,1);
ib_s_all(:,1) = Ib_s * shock_s + Ib_ib_true * nu_all(:,1);
error_s_all(:,1) = A_y * y_s_all(:,1) + A_pi * pi_s_all(:,1) + A_ib * ib_s_all(:,1);

% construct all subsequent shocks and IRFs

for t = 2:T-h+1

    nu_all(1:h,t) = -(A_pi(t:h+t-1,:) * [zeros(t-1,h);Pi_ib_true(1:T-t+1,1:h)] + A_y(t:h+t-1,:) * [zeros(t-1,h);Y_ib_true(1:T-t+1,1:h)] ...
        + A_ib(t:h+t-1,:) * [zeros(t-1,h);Ib_ib_true(1:T-t+1,1:h)])^(-1) * (A_pi(t:h+t-1,:) * pi_s_all(:,t-1) ...
        + A_y(t:h+t-1,:) * y_s_all(:,t-1) + A_ib(t:h+t-1,:) * ib_s_all(:,t-1));
    
    y_s_all(:,t)  = y_s_all(:,t-1) + [zeros(t-1,T);Y_ib_true(1:T-t+1,:)] * nu_all(:,t);
    pi_s_all(:,t) = pi_s_all(:,t-1) + [zeros(t-1,T);Pi_ib_true(1:T-t+1,:)] * nu_all(:,t);
    ib_s_all(:,t) = ib_s_all(:,t-1) + [zeros(t-1,T);Ib_ib_true(1:T-t+1,:)] * nu_all(:,t);
    error_s_all(:,t) = A_y * y_s_all(:,t) + A_pi * pi_s_all(:,t) + A_ib * ib_s_all(:,t);

end

% collect the final ex-post IRFs

y_s_pred  = y_s_all(:,T-h+1);
pi_s_pred = pi_s_all(:,T-h+1);
ib_s_pred = ib_s_all(:,T-h+1);

% compute interest rate surprises

ib_s_surprise = NaN(T-h+1,1);
ib_s_surprise(1) = sum(sqrt(ib_s_all(:,1).^2));
for t = 2:T-h+1
    ib_s_surprise(t) = sum(sqrt((ib_s_all(:,t)-ib_s_all(:,t-1)).^2));
end
ib_s_surprise = [ib_s_surprise;zeros(h-1,1)];
ib_s_surprise = ib_s_surprise ./ ib_s_surprise(1);

end