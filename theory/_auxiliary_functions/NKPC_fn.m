function pi_seq = NKPC_fn(y_seq,shock_s);

global T kappa r_b_SS

pi_seq = zeros(T,1);

pi_seq(T) = kappa * y_seq(T) + shock_s(T);
for t = T-1:-1:1
    pi_seq(t) = kappa * y_seq(t) + 1/(1+r_b_SS) * pi_seq(t+1) + shock_s(t);
end