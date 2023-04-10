function [y_z_optol,pi_z_optpol,ib_z_optpol] = optpol_fn(W_pi,W_y,...
    Theta_pi,Theta_y,Theta_ib,pi_z,y_z,ib_z);

% get optimal shock sequence

nu_z_optpol = - (Theta_pi' * W_pi * Theta_pi + Theta_y' * W_y * Theta_y)^(-1) ...
    * (Theta_pi' * W_pi * pi_z + Theta_y' * W_y * y_z);

% get outcomes

y_z_optol  = y_z + Theta_y * nu_z_optpol;
pi_z_optpol = pi_z + Theta_pi * nu_z_optpol;
ib_z_optpol = ib_z + Theta_ib * nu_z_optpol;

end