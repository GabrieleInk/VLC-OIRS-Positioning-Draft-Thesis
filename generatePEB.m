function crlb = generatePEB(x_scan,y_scan,LED, RIS, N, Phi_FoV, a, rho, Psi, A_pd, T_of,R_pd,p,alp,bet,K0,KN,m,alpha,beta,w,h) % restituisce il limite inferiore per l'errore di precisione e quindi l'errore minimo che si commette 

crlb=zeros(length(x_scan),length(y_scan));

for i=1:length(x_scan)
    parfor j=1:length(y_scan)
        test=[x_scan(i),y_scan(j),0];
        covering_ris=0;
        for cov_index=1:size(RIS,1)
            if isCovered(LED', RIS(cov_index,:)', test', alpha(cov_index), beta(cov_index),N(cov_index,:),w(cov_index),h(cov_index))
                covering_ris=covering_ris+1;
            end
        end
        if covering_ris>1
            crlb_matrix=CRLB_Position(RIS,test,LED,alp,bet,K0,KN,m,R_pd,p,LoS_Contribution(LED', test', Phi_FoV, a, Psi, A_pd, T_of),Phi_FoV,a,N,rho, Psi, A_pd, T_of);
            crlb(i,j)=sqrt(crlb_matrix(1,1)+crlb_matrix(2,2));
        else
            crlb(i,j)=NaN;
        end

    end
end
end