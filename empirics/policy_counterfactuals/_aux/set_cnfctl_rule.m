% NGDP targeting

if cnfctl_ngdp == 1

A_pi = eye(T);
A_y  = eye(T);
for i = 2:T
    A_y(i,i-1) = -1;
end
A_ib = zeros(T,T);

end

% stabilize output

if cnfctl_0y == 1

A_pi = zeros(T,T);
A_y  = eye(T);
A_ib = zeros(T,T);

end

% stabilize interest rates

if cnfctl_0ib == 1

A_pi = zeros(T,T);
A_y  = zeros(T,T);
A_ib = eye(T);

end

% Taylor rule

if cnfctl_tylr == 1

rho_ib = 0.5;

A_pi = (1-rho_ib) * (-1.5 * eye(T));
A_y  = (1-rho_ib) * (-1 * eye(T));
A_ib = eye(T);
for t = 2:T
    A_ib(t,t-1) = -rho_ib;
end

end