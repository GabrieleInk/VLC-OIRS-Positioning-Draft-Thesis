function [Sx,Sz] = Sx_Sz...
    (LED, RIS, PD, alpha, beta,N,w,h)   
% PARAMETERS
% LED =
% RIS =
% PD =
% alpha = Rotation around x [rad]
% beta = Rotation around z [rad]
% N = 
% w = 
% h =
% RETURNS
% h_0 =


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


end