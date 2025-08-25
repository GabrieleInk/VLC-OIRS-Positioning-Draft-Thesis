function plotScenario(LED, RISs, alphas, betas, Ns, ws, hs)

scatter3(LED(1),LED(2),LED(3), 'MarkerFaceColor',		[1 0 0])
hold on

for i=1:size(RISs,1)

if abs(Ns(i,2))==1
rect=[RISs(i,:)+[ws(1)/2 0 hs(1)/2];RISs(i,:)+[-ws(1)/2 0 hs(1)/2];RISs(i,:)+[-ws(1)/2 0 -hs(1)/2];RISs(i,:)+[ws(1)/2 0 -hs(1)/2]];
else
rect=[RISs(i,:)+[ 0 ws(1)/2 hs(1)/2];RISs(i,:)+[ 0 -ws(1)/2 hs(1)/2];RISs(i,:)+[ 0 -ws(1)/2 -hs(1)/2];RISs(i,:)+[ 0 ws(1)/2 -hs(1)/2]];
end

rect1=rect-RISs(i,:);

if Ns(i,1)==1 && Ns(i,2)==0 && Ns(i,3)==0
    theta_y = pi;  
    theta_z = (3/2) * pi;  
elseif Ns(i,1)==-1 && Ns(i,2)==0 && Ns(i,3)==0
    theta_y = pi;  
    theta_z = pi/2;  
elseif Ns(i,1)==0 && Ns(i,2)==1 && Ns(i,3)==0
    theta_y = pi;  
    theta_z = 0;  
elseif Ns(i,1)==0 && Ns(i,2)==-1 && Ns(i,3)==0
    theta_y = pi;  
    theta_z = pi;  
end
  
Ry = [cos(theta_y) 0 -sin(theta_y);
      0 1 0;
      sin(theta_y) 0 cos(theta_y)];

Rz = [cos(theta_z) sin(theta_z) 0;
      -sin(theta_z) cos(theta_z) 0;
      0 0 1];

h_tilt = [1 0 0;
      0 cos(-betas(i)) -sin(-betas(i));
      0 sin(-betas(i)) cos(-betas(i))];

v_tilt = [cos(-alphas(i)) sin(-alphas(i)) 0;
      -sin(-alphas(i)) cos(-alphas(i)) 0;
      0 0 1];

rect2(1,:)=Ry * Rz * rect1(1,:)';
rect2(2,:)=Ry * Rz * rect1(2,:)';
rect2(3,:)=Ry * Rz * rect1(3,:)';
rect2(4,:)=Ry * Rz * rect1(4,:)';

rect3(1,:)=h_tilt * v_tilt * rect2(1,:)';
rect3(2,:)=h_tilt * v_tilt * rect2(2,:)';
rect3(3,:)=h_tilt * v_tilt * rect2(3,:)';
rect3(4,:)=h_tilt * v_tilt * rect2(4,:)';

Ry_back = [cos(-theta_y) 0 -sin(-theta_y);
      0 1 0;
      sin(-theta_y) 0 cos(-theta_y)];

Rz_back = [cos(-theta_z) sin(-theta_z) 0;
      -sin(-theta_z) cos(-theta_z) 0;
      0 0 1];

rect4(1,:)=Rz_back * Ry_back * rect3(1,:)';
rect4(2,:)=Rz_back * Ry_back * rect3(2,:)';
rect4(3,:)=Rz_back * Ry_back * rect3(3,:)';
rect4(4,:)=Rz_back * Ry_back * rect3(4,:)';

rect4=rect4+RISs(i,:);

fill3(rect4(:,1),rect4(:,2),rect4(:,3),[0 0.4470 0.7410])
hold on
end


end