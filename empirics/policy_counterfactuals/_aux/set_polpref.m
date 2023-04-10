lambda_pi_avg = 0.6;  % weight on average inflation deviation
lambda_pi_t   = 1 - lambda_pi_avg;  % weight on current inflation deviation
lambda_y      = 1;    % weight on output deviation

R = 20;      % periods over which we average inflation
delta = 0.1; % smoothing parameter

beta = 1/1.01; % discount factor

% auxiliary matrix

lagweights = exp(-delta * (0:R-1)');
lagweights = lagweights/sum(lagweights);
lagweights_flip = flipud(lagweights);

Pi_bar = zeros(T,T);
for t=1:T
    for s = 1:T
        if s <= t && s > t-R
            Pi_bar(t,s) = lagweights(t+1-s);
        end
    end
end

% construct the weighting matrices

W = zeros(T,T);

for t = 1:T 
    W(t,t) = beta^(t-1);
end

W_pi = lambda_pi_avg * Pi_bar' * W * Pi_bar + lambda_pi_t * W;
W_y  = lambda_y * W;