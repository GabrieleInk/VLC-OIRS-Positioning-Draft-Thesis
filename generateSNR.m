function rem = generateSNR(x_scan, y_scan, LED, RIS, alpha, beta, N, w, h, Phi_FoV, a, rho, Psi, A_pd, T_of,R_pd,p,alp,bet)

rem=zeros(length(x_scan'),length(y_scan'));

for i=1:length(x_scan)
    parfor j=1:length(y_scan)
        test=[x_scan(i),y_scan(j),0];
        signal = 0;
        h0 = LoS_Contribution(LED', test', Phi_FoV, a, Psi, A_pd, T_of);
        signal=signal + R_pd*p*h0;
        for ris_index=1:size(RIS,1)
            [hn,Sx,Sz] = NLoS_Contribution(LED', RIS(ris_index,:)',test', alpha(ris_index), beta(ris_index), Phi_FoV, a, rho, Psi, A_pd, T_of, N(ris_index,:));
            if abs(Sx)<=w(ris_index)/2 && abs(Sz)<=h(ris_index)/2
                signal = signal + R_pd*p*hn;
            end
        end
        rem(i,j)=rem(i,j) + 10*log10((signal^2)/(alp + bet*signal));
    end
end

end