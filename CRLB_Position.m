function crlb = CRLB_Position(RIS,PD,LED,a,b,K0,KN,m,R_pd,p,h0,Phi_FoV,refr_a,Norm,rho, Psi, A_pd, T_of)
d0_real=calculateDistance(LED,PD);
xi=(R_pd*p*h0)*(d0_real^(3 + m));

J=zeros(size(RIS,1)+1,size(RIS,1)+1);
T=zeros(2,size(RIS,1)+1);
for i=1:(size(RIS,1)+1)
    for j=1:i
        if(i==1 && j==1)
        %Elemento in posizione (0,0)
            J(i,j)=J(i,j)+((d0_real^(-5-m))*K0*((3+m)^2)*(xi^2)*((2*a + b^2)*(d0_real^(3+m)) + 2*b*xi))/...
                (2*(a*(d0_real^(3 + m)) + b*xi)^2);
            for k=1:size(RIS,1)
                dn_real=calculateDistance(PD,RIS(k,:));
                s_n=calculateDistance(LED,RIS(k,:));
                [alpha_real,beta_real] = computeTilt(RIS(k,:),LED,PD,Norm(k,:));
                [hn,Sx,Sz] = NLoS_Contribution(LED', RIS(k,:)', PD', alpha_real, beta_real, Phi_FoV, refr_a, rho, Psi, A_pd, T_of, Norm(k,:));
                kappa=(R_pd*p*hn*(dn_real*(dn_real + s_n)^2));
                J(i,j)=J(i,j) + ...
                    (((b^2)*(d0_real^(-8-(2*m)))*KN*((-3-m)^2)*(xi^2))/...
                    (2*(a + b*(kappa/((dn_real*((dn_real + s_n)^2)))+(d0_real^(-3-m))*xi))^2)) + ...
                    ((KN*(6*(3*(d0_real^(-8-(2*m))) + (d0_real^(-8-(2*m)))*m)*(xi^2) + 2*m*(3*(d0_real^(-8-(2*m))) + (d0_real^(-8-(2*m)))*m)*(xi^2)))/...
                    (2*(a + b*(kappa/(dn_real*((dn_real + s_n)^2)) + (d0_real^(-3-m))*xi))));
            end
        elseif(i==j && j~=1)
            %Elemento in posizione (n,n)
            dn_real=calculateDistance(PD,RIS(j-1,:));
            s_n=calculateDistance(LED,RIS(j-1,:));
            [alpha_real,beta_real] = computeTilt(RIS(j-1,:),LED,PD,Norm(j-1,:));
            [hn,Sx,Sz] = NLoS_Contribution(LED', RIS(j-1,:)', PD', alpha_real, beta_real, Phi_FoV, refr_a, rho, Psi, A_pd, T_of, Norm(j-1,:));
            kappa=((R_pd*p*hn)*(dn_real*(dn_real + s_n)^2));
            J(i,j) = J(i,j) + ...
                (((b^2)*KN*((-((2*kappa)/(dn_real*((dn_real+s_n)^3))) - (kappa/((dn_real^2)*((dn_real+s_n)^2))))^2))/...
                (2*(a + b*((kappa/(dn_real*((dn_real+s_n)^2))) + (d0_real^(-3-m))*xi))^2)) + ...
                ((KN*((8*(kappa^2))/((dn_real^2)*(dn_real+s_n)^6) + (8*(kappa^2))/((dn_real^3)*((dn_real+s_n)^5)) + (2*(kappa^2))/((dn_real^4)*((dn_real+s_n)^4))))/...
                (2*(a + b*(kappa/(dn_real*((dn_real+s_n)^2)) + (d0_real^(-3-m))*xi))));
        elseif(j==1 && i~=1)
            dn_real=calculateDistance(PD,RIS(i-1,:));
            s_n=calculateDistance(LED,RIS(i-1,:));
            [alpha_real,beta_real] = computeTilt(RIS(i-1,:),LED,PD,Norm(i-1,:));
            [hn,Sx,Sz] = NLoS_Contribution(LED', RIS(i-1,:)', PD', alpha_real, beta_real, Phi_FoV, refr_a, rho, Psi, A_pd, T_of, Norm(i-1,:));
            kappa=((R_pd*p*hn)*(dn_real*(dn_real + s_n)^2));
            J(i,j) = J(i,j) + ...
                (((b^2)*(d0_real^(-4-m))*KN*(-3-m)*(-((2*kappa)/(dn_real*((dn_real+s_n)^3))) - kappa/((dn_real^2)*((dn_real+s_n)^2)))*xi)/...
                (2*(a+b*(kappa/(dn_real*(dn_real+s_n)^2)+(d0_real^(-3-m))*xi))^2)) + ...
                ((KN*(-((4*(-3*(d0_real^(-4-m))-(d0_real^(-4-m))*m)*kappa*xi)/(dn_real*((dn_real+s_n)^3))) - (2*(-3*(d0_real^(-4-m)) - (d0_real^(-4-m))*m)*kappa*xi)/((dn_real^2)*(dn_real+s_n)^2)))/...
                (2*(a+b*(kappa/(dn_real*(dn_real+s_n)^2)+(d0_real^(-3-m))*xi))));
            J(j,i)=J(i,j);
        end
    end
end
T(1,1)= (PD(1)-LED(1))/(sqrt((PD(1)-LED(1))^2 + (PD(2)-LED(2))^2 + LED(3)^2));
T(2,1)= (PD(2)-LED(2))/(sqrt((PD(1)-LED(1))^2 + (PD(2)-LED(2))^2 + LED(3)^2));
for k=1:size(RIS,1)
    T(1,1+k)= (PD(1)-RIS(k,1))/(sqrt((PD(1)-RIS(k,1))^2 + (PD(2)-RIS(k,2))^2 + RIS(k,3)^2));
    T(2,1+k)= (PD(2)-RIS(k,2))/(sqrt((PD(1)-RIS(k,1))^2 + (PD(2)-RIS(k,2))^2 + RIS(k,3)^2));
end
fim=T*J*T';
crlb=inv(fim);
end