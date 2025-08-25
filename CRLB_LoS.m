function crlb = CRLB_LoS(R_pd,p,h0,d0_real,m,alpha,beta,K)
xi=(R_pd*p*h0)*(d0_real^(3 + m)); %in xi ho la dipendenza dalla potenza e quindi dall'SNR 
crlb=(2*(alpha*(d0_real^(3 + m)) + beta*xi)^2)/...
        ((d0_real^(-5 - m))*K*((3 + m)^2)*(xi^2)*((2*alpha + (beta^2))*(d0_real^(3 + m)) + 2*beta*xi));
end