%% RMSE confrontato con DEB degli stimatori LoS ML e RML per valori diversi di campioni di luce non coerente K, in funzione del SNR

clear variables;
tic
% Random seed
rng(1);

% Load workspace from config.mat
load("config.mat");

% Set entities position

LED     = [x_max/2      , y_max/2      , z_max    ];

SNR_db=10:1:25; % è un vettore di elementi che partono da 10 ed incrementano di 1 fino ad arrivare a 25 (16 elementi)

test_samples=[1,2,3,4,5];
    
RMSE_RML=zeros(length(test_samples),length(SNR_db));
RMSE_ML=zeros(length(test_samples),length(SNR_db));
CRLB=zeros(length(test_samples),length(SNR_db));

times=10000;

d0_r_v=zeros(length(SNR_db),times);
d0_e_v=zeros(length(SNR_db),times);
dist_r=zeros(length(test_samples),length(SNR_db));
dist_e=zeros(length(test_samples),length(SNR_db));

PD_test = [x_max/2 , y_max/2  , 0];
d0_real=calculateDistance(LED,PD_test);

for sample_index=1:length(test_samples)

    h0 = LoS_Contribution(LED', PD_test', Phi_FoV, a, Psi, A_pd, T_of);
    [alpha, beta] = noiseParameters(q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    for test_snr=1:length(SNR_db)
        s=10^(SNR_db(test_snr)/10);
        p=(s*beta + h0*R_pd*sqrt((s*(4*alpha + s*(beta^2)))/((h0^2)*(R_pd^2))))/(2*h0*R_pd);
        
        sum_RML=0;
        sum_ML=0;
        parfor i=1:times % le iterazioni vengono distribuite su più core

            [noises, thermal_var] = noisesEstimation(R_pd*p*h0, q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B, test_samples(sample_index));
            
            mu_0=R_pd*p*h0 + noises;
            mu_0_mean = (sum(mu_0))/test_samples(sample_index);

            d0_r = RML_LoS(R_pd, p, A_pd, G*T_of, m, LED(3), mu_0_mean);
            d0_e = ML_LoS(test_samples(sample_index), thermal_var, q_0, B, I_bg, R_pd, p, A_pd, G*T_of, m, LED(3), mu_0);

            sum_RML=sum_RML + (d0_real - d0_r)^2;
            sum_ML=sum_ML + (d0_real - d0_e)^2;
            
            d0_r_v(test_snr,i)= abs(d0_real - d0_r); % errore assoluto per RML
            d0_e_v(test_snr,i)= abs(d0_real - d0_e); % errore assoluto per ML
        end
        RMSE_RML(sample_index,test_snr)= sqrt(sum_RML/times);
        RMSE_ML(sample_index,test_snr)= sqrt(sum_ML/times);
        CRLB(sample_index,test_snr)=sqrt(CRLB_LoS(R_pd,p,h0,d0_real,m,alpha,beta,test_samples(sample_index)));
    end
    dist_e(sample_index,:)=mean(d0_e_v'); % media dell'errore assoluto per RML
    dist_r(sample_index,:)=mean(d0_r_v'); % media dell'errore assoluto per ML
end

save("figure2.mat","CRLB","RMSE_ML","RMSE_RML","SNR_db","dist_r","dist_e");
elapsedtime2=toc
if isfile("savings.mat")
    load("savings.mat");  % carica la struct "savings"
else
    savings = struct();   % crea la struct se non esiste
end

savings.Figure2_generate = elapsedtime2;  % aggiungi/aggiorna campo

save("savings.mat","savings");