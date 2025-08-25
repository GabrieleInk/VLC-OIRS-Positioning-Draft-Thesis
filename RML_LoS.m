function d_0 = RML_LoS(R, p, A, G, m, q_z, mu_0)
% PARAMETERS
% R = 
% p = 
% A = 
% G = 
% m = 
% q_z = 
% mu_0 = 
% RETURNS
% d_0 = 

xi = (R*p*A*G*(m+1)*(q_z^(m+1)))/(2*pi);
d_0 =(xi/mu_0)^(1/(m+3));
end