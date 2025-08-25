function impulseResponse = estimated_LoS_Contribution...
    (LED, d_SD, Phi_FoV, a, Psi, A_pd, T_of)   
% PARAMETERS
% LED =
% PD =
% Phi_FoV = 
% a =
% Psi = 
% A_pd = 
% T_of = 
% RETURNS
% h_0 =

% Compute lambertian index
m = -(log(2)/log(cosd(Psi))); 
G = (a^2)/((sind(Phi_FoV)^2));
    
impulseResponse = ((m+1)*A_pd * T_of * G * LED(3)^(m+1))/ (2* pi * (d_SD)^(m+3)) ;

end

