clear variables;

% Random seed
rng(1);

% Load workspace from config.mat
load("config.mat");

% Set entities position
LED     = [x_max/2      , y_max/2      , z_max    ];

RIS(1,:)    = [x_max/2      , 0            , z_max/2];
Norm(1,:)=[0 1 0];
alpha(1) = 0;
beta(1) = 0;
w(1)=1;
h(1)=1;

SNR_db=30:1:60;

test_samples=[1,3,5];

RMSE_RML=zeros(length(test_samples),length(SNR_db));
RMSE_ML=zeros(length(test_samples),length(SNR_db));
CRLB=zeros(length(test_samples),length(SNR_db));

times=10000;

d1_r_v=zeros(length(SNR_db),times);
d1_e_v=zeros(length(SNR_db),times);
dist_r=zeros(length(test_samples),length(SNR_db));
dist_e=zeros(length(test_samples),length(SNR_db));

PD_test = [3 , 3 , 0 ];

[alpha_real,beta_real] = computeTilt(RIS(1,:),LED,PD_test,Norm(1,:));
[hn_tilt,Sx_tilt,Sz_tilt] = NLoS_Contribution(LED', RIS(1,:)', PD_test', alpha_real, beta_real, Phi_FoV, a, rho, Psi, A_pd, T_of, Norm(1,:));
[hn,Sx,Sz] = NLoS_Contribution(LED', RIS(1,:)', PD_test', alpha(1), beta(1), Phi_FoV, a, rho, Psi, A_pd, T_of, Norm(1,:));
h0 = LoS_Contribution(LED', PD_test', Phi_FoV, a, Psi, A_pd, T_of);
[alp, bet] = noiseParameters(q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
d1_real=calculateDistance(PD_test,RIS(1,:));
d0_real = calculateDistance(PD_test,LED);
r_1=calculateDistance(LED,RIS(1,:));
for sample_index=1:length(test_samples)
    for test_snr=1:length(SNR_db)
        s=10^(SNR_db(test_snr)/10);
        p=(s*bet + (h0+hn)*R_pd*sqrt((s*(4*alp + s*(bet^2)))/(((h0+hn)^2)*(R_pd^2))))/(2*(h0+hn)*R_pd);
        sum_RML=0;
        sum_ML=0;
        parfor i=1:times
            [noises, thermal_var] = noisesEstimation_v2(R_pd*p*h0 + R_pd*p*hn,R_pd*p*h0, q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B, test_samples(sample_index));
            mu_1 = R_pd*p*hn + noises;
            mu_1_mean = (sum(mu_1))/test_samples(sample_index);

            d1_r = RML_NLoS(rho, R_pd, p, A_pd, G*T_of, m, LED(3), RIS(1,3), r_1, mu_1_mean);
            d1_e = ML_NLoS(rho, test_samples(sample_index), thermal_var, q_0, B, I_bg, R_pd, p, A_pd, G*T_of, m, LED(3), RIS(1,3), r_1, mu_1,d0_real);

            sum_RML=sum_RML + (d1_real-d1_r)^2;
            sum_ML=sum_ML + (d1_real-d1_e)^2;

            d1_r_v(test_snr,i)= abs(d1_real - d1_r);
            d1_e_v(test_snr,i)= abs(d1_real - d1_e);
        end

        RMSE_RML(sample_index,test_snr)= sqrt(sum_RML/times);
        RMSE_ML(sample_index,test_snr)= sqrt(sum_ML/times);
        CRLB(sample_index,test_snr)=sqrt(CRLB_NLoS(R_pd,p,h0,hn_tilt,d0_real,d1_real,r_1,m,alp,bet,test_samples(sample_index)));
    end
    dist_e(sample_index,:)=mean(d1_e_v');
    dist_r(sample_index,:)=mean(d1_r_v');
end


save("figure1_7.mat","CRLB","RMSE_ML","RMSE_RML","SNR_db","dist_r","dist_e");


