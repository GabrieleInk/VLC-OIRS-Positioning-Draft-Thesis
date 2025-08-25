function [impulseResponse,Sx,Sz] = NLoS_Contribution...
    (LED, RIS, PD, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, N)   
% PARAMETERS
% LED =
% RIS =
% PD =
% alpha = Rotation around x [rad]
% beta = Rotation around z [rad]
% Phi_FoV = 
% a =
% Psi = 
% A_pd = 
% T_of = 
% N = Normal vector of the surface
% RETURNS
% h_0 =

% Compute lambertian index
m = -(log(2)/log(cosd(Psi))); 

% (1) Move the reference system in the center of the RIS
PD1=PD-RIS;
LED1=LED-RIS;
RIS1=RIS-RIS;

% Define the angles (in radians)
if N(1)==1 && N(2)==0 && N(3)==0
    theta_y = pi;  
    theta_z = (3/2) * pi;  
elseif N(1)==-1 && N(2)==0 && N(3)==0
    theta_y = pi;  
    theta_z = pi/2;  
elseif N(1)==0 && N(2)==1 && N(3)==0
    theta_y = pi;  
    theta_z = 0;  
elseif N(1)==0 && N(2)==-1 && N(3)==0
    theta_y = pi;  
    theta_z = pi;  
end

% Rotation matrices
  
Ry = [cos(theta_y) 0 -sin(theta_y);
      0 1 0;
      sin(theta_y) 0 cos(theta_y)];

Rz = [cos(theta_z) sin(theta_z) 0;
      -sin(theta_z) cos(theta_z) 0;
      0 0 1];

h_tilt = [1 0 0;
      0 cos(beta) -sin(beta);
      0 sin(beta) cos(beta)];

v_tilt = [cos(alpha) sin(alpha) 0;
      -sin(alpha) cos(alpha) 0;
      0 0 1];

% (2) Apply rotation to have y // N and and z pointing down

PD2=Ry * Rz * PD1;
LED2=Ry * Rz * LED1;
RIS2=Ry * Rz * RIS1;

% (3) Add tilt

PD3=h_tilt * v_tilt * PD2;
LED3=h_tilt * v_tilt * LED2;
RIS3=h_tilt * v_tilt * RIS2;


Sx = ((LED3(2)*PD3(1)) + (LED3(1)*PD3(2)))/(LED3(2) + PD3(2));
Sz = ((LED3(3)*PD3(2)) + (LED3(2)*PD3(3)))/(LED3(2) + PD3(2));

RIS3(1) = RIS3(1) + Sx;     %[Sx, 0, Sz]
RIS3(3) = RIS3(3) + Sz;

h_tilt_back = [1 0 0;
      0 cos(-beta) -sin(-beta);
      0 sin(-beta) cos(-beta)];

v_tilt_back = [cos(-alpha) sin(-alpha) 0;
      -sin(-alpha) cos(-alpha) 0;
      0 0 1];

PD4=v_tilt_back * h_tilt_back * PD3;
LED4=v_tilt_back * h_tilt_back * LED3;
RIS4=v_tilt_back * h_tilt_back * RIS3;

theta_RS = calculateAngle(LED4, RIS4); % angle of irradiance % LED to RIS
Phi_RD = calculateAngle(RIS4, PD4); % angle of incidence % RIS to PD
    
d_SR = calculateDistance(LED4, RIS4); % distance LED to RIS
d_RD = calculateDistance(RIS4, PD4); % distance RIS to PD
    
    % optical concentrator gain
    if Phi_RD >= 0 && Phi_RD <= Phi_FoV
        G = (a^2)/((sind(Phi_FoV)^2));
    else
        G = 0;
    end
    
    impulseResponse = ((rho*(m+1)*A_pd * (cosd(theta_RS)^m) * cosd(Phi_RD) * T_of * G)/ (2* pi * (d_RD+d_SR)^2)) ;

end