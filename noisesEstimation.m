function [noise,t] = noisesEstimation(P_curr, q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B, samples)

var_shoot_first_term = 2*q_0*P_curr*B; % A * s * A * (1/s) % sigma s^2 primo addendo

var_shoot_second_term = 2*q_0*I_bg*B; % caso I_DC = 5 pA % sigma s^2 secondo addendo

var_shoot = var_shoot_first_term + var_shoot_second_term; % sigma s^2

shoot=zeros(samples,1);

% thermal noise
var_thermal = ((8*pi*k_B*T_k*eta*A_pd*I_2*B^2)/G_0 + ((4*pi)^2*k_B*T_k*Gamma*eta^2*A_pd^2*I_3*B^3)/g_m); % sigma t^2

thermal=zeros(samples,1);
for i=1:samples
    shoot(i)= randn(1)*sqrt(var_shoot);
    thermal(i)= randn(1)*sqrt(var_thermal);
    %Remove negative outliers
    % while (P_curr + shoot(i) + thermal(i)) <= 0
    %     shoot(i)= randn(1)*sqrt(var_shoot);
    %     thermal(i)= randn(1)*sqrt(var_thermal);
    % end
end

if P_curr > 0 
    noise = shoot + thermal;
    t=var_thermal;
else
    noise = 0;
    t=var_thermal;
end

end
