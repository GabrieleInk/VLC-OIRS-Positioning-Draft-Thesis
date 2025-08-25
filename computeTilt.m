function [alpha,beta] = computeTilt(RIS,LED,PD_est, N)


LED= LED - RIS;
PD_est=PD_est - RIS;
RIS= RIS-RIS;

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

PD_est =(Ry * Rz * PD_est');
LED = (Ry * Rz * LED');
RIS = (Ry * Rz * RIS');



RS=((LED-RIS)/norm(LED-RIS));
RD=((PD_est-RIS)/norm(PD_est-RIS));
N1=((RS+RD)/sqrt(2+(2*RS'*RD)));
e3=[0 0 1]';
e1=[1 0 0]';
beta = -asin(N1'*e3);
% if N(1)==-1 || N(2)==-1
alpha = -asin(N1'*e1/cos(beta));
% else
% alpha = -asin(N1'*e1/cos(beta));
% end
end

