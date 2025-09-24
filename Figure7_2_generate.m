%% Performance di localizzazione basata sulla legenda di colori (freddo (~1mm) - caldo (~4mm) quindi errore commettibile più alto) con RMSE, flusso luminoso fisso a 1000lm  e 2 iterazioni dell'algoritmo di adaptive beam steering

clear variables;
tic
% Random seed
rng(1);

% Load workspace from config.mat
load("config.mat");

K0 = 50;
KN = 100;

% Set entities position
LED     = [x_max/2      , y_max/2      , z_max    ];

trials=10; % trials=1000;
granularity=20; % granularity=1000;
thresh=10e-3; % soglia
iterations=5;

res=x_max/granularity;
x_scan=0:res:x_max;
y_scan=0:res:y_max;
[x,y]=meshgrid(0:res:x_max,0:res:y_max);

RMSE_1=NaN(granularity+1,granularity+1);
RMSE_2=NaN(granularity+1,granularity+1);
RMSE_3=NaN(granularity+1,granularity+1);

for RIS_config=1:3 % for RIS_config=1:3 per eseguire tutti e tre i casi
    clear RIS Norm alpha beta w h;

    % 4 OIRS
    if RIS_config==1
        RIS=zeros(4,3);
        Norm=zeros(4,3);
        initial_alpha=zeros(4,1);
        initial_beta=zeros(4,1);
        w=zeros(4,1);
        h=zeros(4,1);

        RIS(1,:)    = [x_max/2      , 0            , 1.35];
        Norm(1,:)=[0 1 0];
        initial_alpha(1) = 0;
        initial_beta(1) = 0;
        w(1)=1;
        h(1)=1;

        RIS(2,:)    = [0            , y_max/2      , 1.35];
        Norm(2,:)=[1 0 0];
        initial_alpha(2) = 0;
        initial_beta(2) = 0;
        w(2)=1;
        h(2)=1;

        RIS(3,:)    = [x_max        , y_max/2      , 1.35];
        Norm(3,:)=[-1 0 0];
        initial_alpha(3) = 0;
        initial_beta(3) = 0;
        w(3)=1;
        h(3)=1;

        RIS(4,:)    = [x_max/2      , y_max        , 1.35];
        Norm(4,:)=[0 -1 0];
        initial_alpha(4) = 0;
        initial_beta(4) = 0;
        w(4)=1;
        h(4)=1;
    
    % 8 OIRS
    elseif RIS_config==2
        RIS=zeros(8,3);
        Norm=zeros(8,3);
        initial_alpha=zeros(8,1);
        initial_beta=zeros(8,1);
        w=zeros(8,1);
        h=zeros(8,1);

        RIS(1,:)    = [(x_max/5)*3      , 0            , 1.35];
        Norm(1,:)=[0 1 0];
        initial_alpha(1) = 0;
        initial_beta(1) = 0;
        w(1)=1;
        h(1)=1;

        RIS(5,:)    = [(x_max/5)*2      , 0            , 1.35];
        Norm(5,:)=[0 1 0];
        initial_alpha(5) = -0;
        initial_beta(5) = 0;
        w(5)=1;
        h(5)=1;

        RIS(2,:)    = [0            , (y_max/5)*3      , 1.35];
        Norm(2,:)=[1 0 0];
        initial_alpha(2) = -0;
        initial_beta(2) = 0;
        w(2)=1;
        h(2)=1;

        RIS(6,:)    = [0            , (y_max/5)*2      , 1.35];
        Norm(6,:)=[1 0 0];
        initial_alpha(6) = 0;
        initial_beta(6) = 0;
        w(6)=1;
        h(6)=1;

        RIS(3,:)    = [x_max        , (y_max/5)*3      , 1.35];
        Norm(3,:)=[-1 0 0];
        initial_alpha(3) = 0;
        initial_beta(3) = 0;
        w(3)=1;
        h(3)=1;

        RIS(7,:)    = [x_max        , (y_max/5)*2      , 1.35];
        Norm(7,:)=[-1 0 0];
        initial_alpha(7) = -0;
        initial_beta(7) = 0;
        w(7)=1;
        h(7)=1;

        RIS(4,:)    = [(x_max/5)*3      , y_max        , 1.35];
        Norm(4,:)=[0 -1 0];
        initial_alpha(4) = -0;
        initial_beta(4) = 0;
        w(4)=1;
        h(4)=1;

        RIS(8,:)    = [(x_max/5)*2      , y_max        , 1.35];
        Norm(8,:)=[0 -1 0];
        initial_alpha(8) = 0;
        initial_beta(8) = 0;
        w(8)=1;
        h(8)=1;
    
    % 12 OIRS   
    elseif RIS_config==3
        RIS=zeros(12,3);
        Norm=zeros(12,3);
        initial_alpha=zeros(12,1);
        initial_beta=zeros(12,1);
        w=zeros(12,1);
        h=zeros(12,1);

        RIS(1,:)    = [(x_max/10)*7      , 0            , 1.35];
        Norm(1,:)=[0 1 0];
        initial_alpha(1) = 0;
        initial_beta(1) = 0;
        w(1)=1;
        h(1)=1;

        RIS(5,:)    = [(x_max/10)*5      , 0            , 1.35];
        Norm(5,:)=[0 1 0];
        initial_alpha(5) = -0;
        initial_beta(5) = 0;
        w(5)=1;
        h(5)=1;

        RIS(9,:)    = [(x_max/10)*3      , 0            , 1.35];
        Norm(9,:)=[0 1 0];
        initial_alpha(9) = -0;
        initial_beta(9) = 0;
        w(9)=1;
        h(9)=1;

        RIS(2,:)    = [0            , (y_max/10)*7      , 1.35];
        Norm(2,:)=[1 0 0];
        initial_alpha(2) = -0;
        initial_beta(2) = 0;
        w(2)=1;
        h(2)=1;

        RIS(6,:)    = [0            , (y_max/10)*5      , 1.35];
        Norm(6,:)=[1 0 0];
        initial_alpha(6) = 0;
        initial_beta(6) = 0;
        w(6)=1;
        h(6)=1;

        RIS(10,:)    = [0            , (y_max/10)*3      , 1.35];
        Norm(10,:)=[1 0 0];
        initial_alpha(10) = 0;
        initial_beta(10) = 0;
        w(10)=1;
        h(10)=1;

        RIS(3,:)    = [x_max        , (y_max/10)*7      , 1.35];
        Norm(3,:)=[-1 0 0];
        initial_alpha(3) = 0;
        initial_beta(3) = 0;
        w(3)=1;
        h(3)=1;

        RIS(7,:)    = [x_max        , (y_max/10)*5      , 1.35];
        Norm(7,:)=[-1 0 0];
        initial_alpha(7) = -0;
        initial_beta(7) = 0;
        w(7)=1;
        h(7)=1;

        RIS(11,:)    = [x_max        , (y_max/10)*3      , 1.35];
        Norm(11,:)=[-1 0 0];
        initial_alpha(11) = -0;
        initial_beta(11) = 0;
        w(11)=1;
        h(11)=1;

        RIS(4,:)    = [(x_max/10)*7      , y_max        , 1.35];
        Norm(4,:)=[0 -1 0];
        initial_alpha(4) = -0;
        initial_beta(4) = 0;
        w(4)=1;
        h(4)=1;

        RIS(8,:)    = [(x_max/10)*5      , y_max        , 1.35];
        Norm(8,:)=[0 -1 0];
        initial_alpha(8) = 0;
        initial_beta(8) = 0;
        w(8)=1;
        h(8)=1;

        RIS(12,:)    = [(x_max/10)*3      , y_max        , 1.35];
        Norm(12,:)=[0 -1 0];
        initial_alpha(12) = 0;
        initial_beta(12) = 0;
        w(12)=1;
        h(12)=1;

    end


    [alp, bet] = noiseParameters(q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    for x_index=1:granularity/2 + 1
        parfor y_index=1:granularity/2 + 1
        %for y_index=1:x_index
            test=[x_scan(x_index),y_scan(y_index),0];
            %%% START MONTECARLO %%%
            result=zeros(length(trials),1);
            covering_ris=0;
            for cov_index=1:size(RIS,1)
                if isCovered(LED', RIS(cov_index,:)', test', initial_alpha(cov_index), initial_beta(cov_index),Norm(cov_index,:),w(cov_index),h(cov_index))
                    covering_ris=covering_ris+1;
                end
            end
            if covering_ris>1
                for i=1:trials

                    %%% LOCALIZATION ALGORITHM WITH RELAXED ESTIMATORS %%%
                    alpha=initial_alpha;
                    beta=initial_beta;

                    PD_est=[0 0 0];


                    anchors=zeros(1,13);
                    h0 = LoS_Contribution(LED', test', Phi_FoV, a, Psi, A_pd, T_of);
                    [noises, thermal_var] = noisesEstimation(R_pd*p*h0, q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B, K0);
                    mu_0 = R_pd*p*h0 + noises;
                    mu_0_mean = (sum(mu_0))/K0;
                    anchors(1,1) = RML_LoS(R_pd, p, A_pd, G*T_of, m, LED(3), mu_0_mean);
                    h0_e=estimated_LoS_Contribution(LED', anchors(1,1), Phi_FoV, a, Psi, A_pd, T_of);
                    for j=1:iterations
                        for k=1:size(RIS,1)
                            [hn,Sx,Sz] = NLoS_Contribution(LED', RIS(k,:)', test', alpha(k), beta(k), Phi_FoV, a, rho, Psi, A_pd, T_of, Norm(k,:));
                            if abs(Sx) > w(k)/2 || abs(Sz) > h(k)/2 || hn == 0
                                anchors(1,1+k) = 0;
                            else
                                [noises, thermal_var] = noisesEstimation(R_pd*p*hn + R_pd*p*h0, q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B, KN);
                                mu_n = R_pd*p*hn + R_pd*p*(h0 - h0_e) + noises;
                                mu_n_mean = (sum(mu_n))/KN;
                                if mu_n_mean>I_bg
                                    r_n=calculateDistance(LED,RIS(k,:));
                                    kappa = ((R_pd*p*rho*(m+1)*A_pd * T_of * G*((LED(3)-RIS(k,3))^m)*RIS(k,3))/ (2* pi));
                                    anchors(1,1+k) = RML_NLoS(rho, R_pd, p, A_pd, G*T_of, m, LED(3), RIS(k,3), r_n, mu_n_mean);
                                else
                                    anchors(1,1+k)=0;
                                end
                            end
                        end
                        distance_vect=[anchors(1,1)];
                        entity_vect=[LED];
                        %d0_real=calculateDistance(LED,test);
                        weight_vect = [1/CRLB_LoS(R_pd,p,h0,anchors(1,1),m,alp,bet,K0)];
                        for l=1:12
                            if anchors(1,1+l)>0
                                distance_vect=[distance_vect,anchors(1,1+l)];
                                entity_vect(length(distance_vect),:) = RIS(l,:);
                                [hn,Sx,Sz] = NLoS_Contribution(LED', RIS(l,:)', test', alpha(l), beta(l), Phi_FoV, a, rho, Psi, A_pd, T_of, Norm(l,:));
                                r_n=calculateDistance(LED,RIS(l,:));
                                %dn_real=calculateDistance(RIS(l,:),test);
                                weight_vect=[weight_vect,1/CRLB_NLoS(R_pd,p,h0,hn,anchors(1,1),anchors(1,1+l),r_n,m,alp,bet,KN)];
                            end
                        end

                        if length(distance_vect)>2
                            if j==1
                                PD_est = WLS(entity_vect, distance_vect, 0, weight_vect);
                            end
                            PD_est = WILS(entity_vect, distance_vect, PD_est(1:2), 0, weight_vect);
                            %Update tilt
                            for k=1:size(RIS,1)
                                [alpha(k),beta(k)] = computeTilt(RIS(k,:),LED,PD_est,Norm(k,:));
                            end
                        end
                    end
                    result(i)= (norm(test-PD_est))^2;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if RIS_config==1
                    RMSE_1(x_index,y_index)=sqrt(mean(result));
                elseif RIS_config==2
                    RMSE_2(x_index,y_index)=sqrt(mean(result));
                else
                    RMSE_3(x_index,y_index)=sqrt(mean(result));
                end
            end
        end
    end
    if RIS_config==1
        save("figure7_2_1.mat","RMSE_1","x","y");
    elseif RIS_config==2
        save("figure7_2_2.mat","RMSE_2","x","y");
    else
        save("figure7_2_3.mat","RMSE_3","x","y");
    end

end
elapsedtime7_2=toc
if isfile("savings.mat")
    load("savings.mat");  % carica la struct "savings"
else
    savings = struct();   % crea la struct se non esiste
end

savings.Figure7_2_generate = elapsedtime7_2;  % aggiungi/aggiorna campo

save("savings.mat","savings");
