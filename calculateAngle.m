
function angle_deg = calculateAngle(u, v)
% PARAMETERS
% u = 
% v = 

d = calculateDistance(u,v);

x = abs(u(3)-v(3)); % d * cos(angle)

angle = acos(x/d);

angle_deg = rad2deg(angle);
end