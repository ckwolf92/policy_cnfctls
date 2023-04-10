function [tau_seq,b_g_seq] = gov_bc_fn(ib_seq,pi_seq,y_seq);

global T r_b_SS y_tax Y_SS B_SS Trans_SS

b_g_seq = zeros(T,1);

tau_seq = 1/Trans_SS * (y_tax * Y_SS * y_seq + B_SS * b_g_seq - (1 + r_b_SS) * B_SS * [0;ib_seq(1:end-1)] ...
    + (1 + r_b_SS) * B_SS * pi_seq - (1 + r_b_SS) * B_SS * [0;b_g_seq(1:end-1)]);

end