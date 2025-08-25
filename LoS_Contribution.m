function impulseResponse = LoS_Contribution...
    (LED, PD, Phi_FoV, a, Psi, A_pd, T_of)   
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

% Compute distance and angle
theta_SD = calculateAngle(LED, PD);
d_SD = calculateDistance(LED, PD);

% Compute optical concentrator gain
if theta_SD >= 0 && theta_SD <= Phi_FoV
    G = (a^2)/((sind(Phi_FoV)^2));
else
    G = 0;
end
    
impulseResponse = (((m+1)*A_pd * ((cosd(theta_SD)^m) * cosd(theta_SD) * T_of * G) )/ (2* pi * (d_SD)^2)) ;

end

