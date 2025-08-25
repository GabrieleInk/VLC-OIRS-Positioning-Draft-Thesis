function crlb = CRLB_NLoS(R_pd,p,h0,hn,d0_real,d1_real,r_1,m,alpha,beta,Kn)
xi=(R_pd*p*h0)*(d0_real^(3 + m)); %in xi ho la dipendenza dalla potenza e quindi dall'SNR  
kappa=(R_pd*p*hn*(d1_real*(d1_real + r_1)^2));
crlb=(2*(d1_real^5)*((d1_real + r_1)^8)*(alpha + (beta*kappa)/(d1_real*((d1_real + r_1)^2)) + (d0_real^(-3 - m))*beta*xi)^2)/...
                   (Kn*((3*d1_real + r_1)^2)*(kappa^2)*(2*d1_real*((d1_real + r_1)^2)*alpha + d1_real*((d1_real + r_1)^2)*(beta^2) + 2*beta*(kappa + d1_real*(d0_real^(-3 - m))*((d1_real + r_1)^2)*xi)));

end