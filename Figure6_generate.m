%% RMSE confrontato con PEB degli stimatori IWLS e ILS per due diversi valori di flusso luminoso p (1000 e 3000 lumens), in funzione delle iterazioni dell'algoritmo di adaptive beam steering

clear variables;
tic
% Random seed
rng(1);

% Load workspace from config.mat
load("config.mat");

% Set entities position
LED     = [x_max/2      , y_max/2      , z_max    ];

trials=10000; % trials=10000;

iterations=5;

K0 = 50;
KN = 100;

RIS=zeros(4,3); % matrice degli specchietti OIRS, sono (x,y) x specchietti e y 3 coordinate
Norm=zeros(4,3);% matrice delle normali degli specchietti OIRS
initial_alpha=zeros(4,1);
initial_beta=zeros(4,1);
alpha=zeros(4,1);
beta=zeros(4,1);
w=zeros(4,1);
h=zeros(4,1);

RIS(1,:)    = [x_max/2      , 0            , 1.5];
Norm(1,:)=[0 1 0];
initial_alpha(1) = 0;
initial_beta(1) = 0;
w(1)=0.6; % lunghezza dell'OIRS
h(1)=0.6; % altezza dell'OIRS

RIS(2,:)    = [0            , y_max/2      , 1.5];
Norm(2,:)=[1 0 0];
initial_alpha(2) = 0;
initial_beta(2) = 0;
w(2)=0.6;
h(2)=0.6;

RIS(3,:)    = [x_max        , y_max/2      , 1.5];
Norm(3,:)=[-1 0 0];
initial_alpha(3) = 0;
initial_beta(3) = 0;
w(3)=0.6;
h(3)=0.6;

RIS(4,:)    = [x_max/2      , y_max        , 1.5];
Norm(4,:)=[0 -1 0];
initial_alpha(4) = 0;
initial_beta(4) = 0;
w(4)=0.6;
h(4)=0.6;

[alp, bet] = noiseParameters(q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
test=[3 3 0];

lumen_level1 = 1000;
p1 = lumen_level1 / 683;      % Transmission power 6000 [Lumens] -> [Watt]

result=zeros(trials,iterations,4);

%%% COMPUTE PEB %%%
crlb_matrix=CRLB_Position(RIS,test,LED,alp,bet,K0,KN,m,R_pd,p1,LoS_Contribution(LED', test', Phi_FoV, a, Psi, A_pd, T_of),Phi_FoV,a,Norm,rho, Psi, A_pd, T_of);
crlb_1=ones(iterations,1)*sqrt(crlb_matrix(1,1)+crlb_matrix(2,2));
%%%%%%%%%%%%%%%%%%%
%[1,"RML ILS"][2,"ML ILS"][3,"RML IWLS"][4,"ML IWLS"]%
for index_met=1:2:3
    %%% START MONTECARLO %%%
    for i=1:trials

        %%% LOCALIZATION ALGORITHM WITH RELAXED ESTIMATORS %%%
        alpha=initial_alpha;
        beta=initial_beta;

        PD_est=[2.5 2.5 0];
        covering_ris=0;
        for cov_index=1:size(RIS,1)
            if isCovered(LED', RIS(cov_index,:)', test', alpha(cov_index), beta(cov_index),Norm(cov_index,:),w(cov_index),h(cov_index))
                covering_ris=covering_ris+1;
            end
        end
        if covering_ris>1
            anchors=zeros(1,5);

            h0 = LoS_Contribution(LED', test', Phi_FoV, a, Psi, A_pd, T_of);
            [noises, thermal_var] = noisesEstimation(R_pd*p1*h0, q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B, K0);
            mu_0 = R_pd*p1*h0 + noises;
            mu_0_mean = (sum(mu_0))/K0;

            if index_met == 1 || index_met == 3
                anchors(1,1) = RML_LoS(R_pd, p1, A_pd, G*T_of, m, LED(3), mu_0_mean);
            else
                anchors(1,1) = ML_LoS(K0, thermal_var, q_0, B, I_bg,R_pd, p1, A_pd, G*T_of, m, LED(3), mu_0);
            end
            h0_e=estimated_LoS_Contribution(LED', anchors(1,1), Phi_FoV, a, Psi, A_pd, T_of);
            for j=1:iterations
                for k=1:size(RIS,1)
                    [hn,Sx,Sz] = NLoS_Contribution(LED', RIS(k,:)', test', alpha(k), beta(k), Phi_FoV, a, rho, Psi, A_pd, T_of, Norm(k,:));
                    if abs(Sx) > w(k)/2 || abs(Sz) > h(k)/2 || hn == 0
                        anchors(1,1+k) = 0;
                    else
                        [noises, thermal_var] = noisesEstimation(R_pd*p1*hn + R_pd*p1*h0, q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B, KN);
                        
                        mu_n = R_pd*p1*hn + R_pd*p1*(h0 - h0_e) + noises;
                        mu_n_mean = (sum(mu_n))/KN;
                        r_n=calculateDistance(LED,RIS(k,:));
                        kappa = ((R_pd*p1*rho*(m+1)*A_pd * T_of * G*((LED(3)-RIS(k,3))^m)*RIS(k,3))/ (2* pi));
                        if index_met == 1 || index_met == 3
                            anchors(1,1+k) = RML_NLoS(rho, R_pd, p1, A_pd, G*T_of, m, LED(3), RIS(k,3), r_n, mu_n_mean);
                        else
                            anchors(1,1+k) = ML_NLoS(rho, KN, thermal_var, q_0, B, I_bg, R_pd, p1, A_pd, G*T_of, m, LED(3), RIS(k,3), r_n, mu_n,anchors(1,1));
                        end

                    end
                end
                distance_vect=[anchors(1,1)];
                entity_vect=[LED];
                if index_met == 3 || index_met == 4
                    %d0_real=calculateDistance(LED,test);
                    weight_vect = [1/CRLB_LoS(R_pd,p1,h0,anchors(1,1),m,alp,bet,K0)];
                end
                for l=1:4
                    if anchors(1,1+l)>0
                        distance_vect=[distance_vect,anchors(1,1+l)];
                        entity_vect(length(distance_vect),:) = RIS(l,:);
                        if index_met == 3 || index_met == 4
                            [hn,Sx,Sz] = NLoS_Contribution(LED', RIS(l,:)', test', alpha(l), beta(l), Phi_FoV, a, rho, Psi, A_pd, T_of, Norm(l,:));
                            r_n=calculateDistance(LED,RIS(l,:));
                            %dn_real=calculateDistance(RIS(l,:),test);
                            weight_vect=[weight_vect,1/CRLB_NLoS(R_pd,p1,h0,hn,anchors(1,1),anchors(1,1+l),r_n,m,alp,bet,KN)];
                        end
                    end
                end

                if length(distance_vect)>2
                    if index_met == 1 || index_met == 2
                        PD_est = ILS(entity_vect, distance_vect, PD_est(1:2), 0);
                    else
                        PD_est = WILS(entity_vect, distance_vect, PD_est(1:2), 0, weight_vect);
                    end
                    
                    %Update tilt
                    for k=1:size(RIS,1)
                        [alpha(k),beta(k)] = computeTilt(RIS(k,:),LED,PD_est,Norm(k,:)); 
                    end
                end
                result(i,j,index_met)= (norm(test-PD_est'))^2;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

rmse_1=sqrt((reshape(sum(result),[5,4])./trials)');

lumen_level2 = 3000;
p2 = lumen_level2 / 683;      % Transmission power 6000 [Lumens] -> [Watt]

result=zeros(trials,iterations,4);

%%% COMPUTE PEB %%%
crlb_matrix=CRLB_Position(RIS,test,LED,alp,bet,K0,KN,m,R_pd,p2,LoS_Contribution(LED', test', Phi_FoV, a, Psi, A_pd, T_of),Phi_FoV,a,Norm,rho, Psi, A_pd, T_of);
crlb_2=ones(iterations,1)*sqrt(crlb_matrix(1,1)+crlb_matrix(2,2));
%%%%%%%%%%%%%%%%%%%
%[1,"RML ILS"][2,"ML ILS"][3,"RML IWLS"][4,"ML IWLS"]%
for index_met=1:2:3
    %%% START MONTECARLO %%%
    for i=1:trials

        %%% LOCALIZATION ALGORITHM WITH RELAXED ESTIMATORS %%%
        alpha=initial_alpha;
        beta=initial_beta;

        PD_est=[2.5 2.5 0];
        covering_ris=0;
        for cov_index=1:size(RIS,1)
            if isCovered(LED', RIS(cov_index,:)', test', alpha(cov_index), beta(cov_index),Norm(cov_index,:),w(cov_index),h(cov_index))
                covering_ris=covering_ris+1;
            end
        end
        if covering_ris>1
            anchors=zeros(1,5);

            h0 = LoS_Contribution(LED', test', Phi_FoV, a, Psi, A_pd, T_of);
            [noises, thermal_var] = noisesEstimation(R_pd*p2*h0, q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B, K0);
            mu_0 = R_pd*p2*h0 + noises;
            mu_0_mean = (sum(mu_0))/K0;

            if index_met == 1 || index_met == 3
                anchors(1,1) = RML_LoS(R_pd, p2, A_pd, G*T_of, m, LED(3), mu_0_mean);
            else
                anchors(1,1) = ML_LoS(K0, thermal_var, q_0, B, I_bg,R_pd, p2, A_pd, G*T_of, m, LED(3), mu_0);
            end
            h0_e=estimated_LoS_Contribution(LED', anchors(1,1), Phi_FoV, a, Psi, A_pd, T_of);
            for j=1:iterations
                for k=1:size(RIS,1)
                    [hn,Sx,Sz] = NLoS_Contribution(LED', RIS(k,:)', test', alpha(k), beta(k), Phi_FoV, a, rho, Psi, A_pd, T_of, Norm(k,:));
                    if abs(Sx) > w(k)/2 || abs(Sz) > h(k)/2 || hn == 0
                        anchors(1,1+k) = 0;
                    else
                        [noises, thermal_var] = noisesEstimation(R_pd*p2*hn + R_pd*p2*h0, q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B, KN);
                        mu_n = R_pd*p2*hn + R_pd*p2*(h0 - h0_e) + noises;
                        mu_n_mean = (sum(mu_n))/KN;
                        r_n=calculateDistance(LED,RIS(k,:));
                        kappa = ((R_pd*p2*rho*(m+1)*A_pd * T_of * G*((LED(3)-RIS(k,3))^m)*RIS(k,3))/ (2* pi));
                        if index_met == 1 || index_met == 3
                            anchors(1,1+k) = RML_NLoS(rho, R_pd, p2, A_pd, G*T_of, m, LED(3), RIS(k,3), r_n, mu_n_mean);
                        else
                            anchors(1,1+k) = ML_NLoS(rho, KN, thermal_var, q_0, B, I_bg, R_pd, p2, A_pd, G*T_of, m, LED(3), RIS(k,3), r_n, mu_n,anchors(1,1));
                        end

                    end
                end
                distance_vect=[anchors(1,1)];
                entity_vect=[LED];
                if index_met == 3 || index_met == 4
                    %d0_real=calculateDistance(LED,test);
                    weight_vect = [1/CRLB_LoS(R_pd,p2,h0,anchors(1,1),m,alp,bet,K0)];
                end
                for l=1:4
                    if anchors(1,1+l)>0
                        distance_vect=[distance_vect,anchors(1,1+l)];
                        entity_vect(length(distance_vect),:) = RIS(l,:);
                        if index_met == 3 || index_met == 4
                            [hn,Sx,Sz] = NLoS_Contribution(LED', RIS(l,:)', test', alpha(l), beta(l), Phi_FoV, a, rho, Psi, A_pd, T_of, Norm(l,:));
                            r_n=calculateDistance(LED,RIS(l,:));
                            %dn_real=calculateDistance(RIS(l,:),test);
                            weight_vect=[weight_vect,1/CRLB_NLoS(R_pd,p2,h0,hn,anchors(1,1),anchors(1,1+l),r_n,m,alp,bet,KN)];
                        end
                    end
                end

                if length(distance_vect)>2
                    if index_met == 1 || index_met == 2
                        PD_est = ILS(entity_vect, distance_vect, PD_est(1:2), 0);
                    else
                        PD_est = WILS(entity_vect, distance_vect, PD_est(1:2), 0, weight_vect);
                    end
                    
                    %Update tilt
                    for k=1:size(RIS,1)
                        [alpha(k),beta(k)] = computeTilt(RIS(k,:),LED,PD_est,Norm(k,:)); 
                    end
                end
                result(i,j,index_met)= (norm(test-PD_est'))^2;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

rmse_2=sqrt((reshape(sum(result),[5,4])./trials)');

save("figure6.mat","rmse_1","rmse_2","crlb_1","crlb_2","iterations");
elapsedtime6=toc
if isfile("savings.mat")
    load("savings.mat");  % carica la struct "savings"
else
    savings = struct();   % crea la struct se non esiste
end

savings.Figure6_generate = elapsedtime6;  % aggiungi/aggiorna campo

save("savings.mat","savings");

