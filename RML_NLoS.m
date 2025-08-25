
function d_n = RML_NLoS(ro, R, p, A, G, m, q_z, wn_z, r_n, mu_n)
% PARAMETERS
% ro = 
% R = 
% p = 
% A = 
% G = 
% m = 
% q_z = 
% wn_z = 
% r_n = 
% mu_n = 
% RETURNS
% d_n = 

kappa = (R*p*ro*A*G*(m+1)*((q_z - wn_z)^m)*wn_z)/(2*pi*(r_n^m));
root_arg = 27*kappa*(mu_n^2) + 3*sqrt(3)*sqrt(kappa*(mu_n^4)*(27*kappa + 4*mu_n*(r_n^3))) + 2*(mu_n^3)*(r_n^3);
gamma = root_arg^(1/3);

d_n = (gamma)/(3*(2^(1/3))*mu_n) + ((2^(1/3))*mu_n*(r_n^2))/(3*gamma) - (2*r_n)/(3);
end