function [alpha,beta] = noiseParameters( q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B)

var_shoot_second_term = 2*q_0*I_bg*B; 

var_thermal = ((8*pi*k_B*T_k*eta*A_pd*I_2*B^2)/G_0 + ((4*pi)^2*k_B*T_k*Gamma*eta^2*A_pd^2*I_3*B^3)/g_m);

alpha = (var_shoot_second_term + var_thermal);
beta= 2*q_0*B;

end
