function [y_z_cnfctl,pi_z_cnfctl,ib_z_cnfctl,nu_z,policy_error] = cnfctl_fn(A_y,A_pi,A_ib,...
    Theta_pi,Theta_y,Theta_ib,pi_z,y_z,ib_z);

% get optimal shock sequence

nu_z = - ((A_pi * Theta_pi + A_y * Theta_y + A_ib * Theta_ib)' * (A_pi * Theta_pi + A_y * Theta_y + A_ib * Theta_ib))^(-1) ...
                * ((A_pi * Theta_pi + A_y * Theta_y + A_ib * Theta_ib)' * (A_pi * pi_z + A_y * y_z + A_ib * ib_z));

% get outcomes

y_z_cnfctl  = y_z + Theta_y * nu_z;
pi_z_cnfctl = pi_z + Theta_pi * nu_z;
ib_z_cnfctl = ib_z + Theta_ib * nu_z;

% get policy error

policy_error = A_y * y_z_cnfctl + A_pi * pi_z_cnfctl + A_ib * ib_z_cnfctl;