function d_0 = ML_LoS(K_0, sigma_T_2, q, B, I_1, R, p, A, G, m, q_z, mu_0)
% PARAMETERS
% K_0 = 
% sigma_T_2 = 
% q =  
% B = 
% I_1 = 
% R = 
% p = 
% A = 
% G = 
% m = 
% q_z = 
% mu_0 = 
% RETURNS
% d_0 = 

alpha= (sigma_T_2) + ((2*q*I_1*B));
beta = (2*q*B); 

%TODO: check condition on mathematica
% if ((alpha/beta) + mu_0) <= sqrt(alpha + (alpha/beta)^2) 
%      error("Condition on exact estimator not observed");
% end

xi = (R*p*A*G*(m+1)*(q_z^(m+1)))/(2*pi);

%With mu_0 negative and positive
%d_0=(2^(-(1)/(m+3)))*(xi*(sqrt((4*(alpha+beta*mu_0)^2+beta^4)/(-alpha*beta+2*alpha*mu_0+beta*(mu_0^2))^2)+...
%    (2*alpha+(beta^2))/(-alpha*beta+2*alpha*mu_0+beta*(mu_0^2))))^((1)/(m+3));
somma=sum(mu_0);
somma2=sum(mu_0.*mu_0);

%g=sqrt(((4*(K_0^2)*(alpha^2) + (beta^4)*(K_0^2) + 8*K_0*alpha*beta*somma + 4*K_0*(beta^2)*somma2)*xi^2)/...
%    ((-(alpha*beta*K_0) + 2*alpha*somma + beta*somma2)^2));

g=sqrt((K_0*(4*(alpha^2)*K_0 + (beta^4)*K_0 + 8*alpha*beta*somma + 4*(beta^2)*somma2)*(xi^2))...
    /(-(alpha*beta*K_0) + 2*alpha*somma + beta*somma2)^2);

%d_0=(-((K_0*2*alpha*xi + (beta^2)*K_0*xi + beta*somma2*g + alpha*(-(beta*K_0) + 2*somma)*g)/...
 %   (2*alpha*beta*K_0 - 4*alpha*somma - 2*beta*somma2)))^(1/(3+m));

d_0=(g/2 - ((2*alpha + beta^2)*K_0*xi)/(2*alpha*beta*K_0 - 4*alpha*somma - 2*beta*somma2))^(1/(3 + m));
end