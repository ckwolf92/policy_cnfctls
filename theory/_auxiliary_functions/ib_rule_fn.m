function ib_seq = ib_rule_fn(pi_seq,y_seq,m_seq);

global T phi_pi phi_y phi_i

ib_seq = NaN(T,1);

ib_seq(1) = (1-phi_i) * (phi_pi * pi_seq(1) + phi_y * y_seq(1) + m_seq(1));
for t = 2:T
    ib_seq(t) = phi_i * ib_seq(t-1) + (1-phi_i) * (phi_pi * pi_seq(t) + phi_y * y_seq(t) + m_seq(t));
end