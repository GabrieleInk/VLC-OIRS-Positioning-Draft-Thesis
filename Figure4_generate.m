clear variables;
tic
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

span=0:0.01:y_max;

test_samples=[5,10,20,40,80];

RMSE_RML=zeros(length(test_samples),length(span));
RMSE_ML=zeros(length(test_samples),length(span));
CRLB=zeros(length(test_samples),length(span));

times=1; % times=10000

d1_r_v=zeros(length(span),times);
d1_e_v=zeros(length(span),times);
dist_r=zeros(length(test_samples),length(span));
dist_e=zeros(length(test_samples),length(span));

PD_test = [x_max/2, 0, 0];

for sample_index=1:length(test_samples)

    [alpha, beta] = noiseParameters(q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    for test_span=1:length(span)
        PD_test = [x_max/2, span(test_span), 0];
        [hn,Sx,Sz] = NLoS_Contribution(LED', RIS(1,:)', PD_test', alpha(1), beta(1), Phi_FoV, a, rho, Psi, A_pd, T_of, Norm(1,:));
        h0 = LoS_Contribution(LED', PD_test', Phi_FoV, a, Psi, A_pd, T_of);
        d1_real=calculateDistance(PD_test,RIS(1,:));
        d0_real = calculateDistance(PD_test,LED);
        r_1=calculateDistance(LED,RIS(1,:));
        sum_RML=0;
        sum_ML=0;
        if isCovered(LED', RIS(1,:)', PD_test', alpha(1), beta(1),Norm(1,:),w(1),h(1))

            parfor i=1:times

                [noises, thermal_var] = noisesEstimation(R_pd*p*h0 + R_pd*p*h0, q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B, test_samples(sample_index));
                mu_1 = R_pd*p*hn + noises;
                mu_1_mean = (sum(mu_1))/test_samples(sample_index);

                d1_r = RML_NLoS(rho, R_pd, p, A_pd, G*T_of, m, LED(3), RIS(1,3), r_1, mu_1_mean);
                d1_e = ML_NLoS(rho, test_samples(sample_index), thermal_var, q_0, B, I_bg, R_pd, p, A_pd, G*T_of, m, LED(3), RIS(1,3), r_1, mu_1,d0_real);

                sum_RML=sum_RML + (d1_real-d1_r)^2;
                sum_ML=sum_ML + (d1_real-d1_e)^2;

                d1_r_v(test_span,i)= abs(d1_real - d1_r);
                d1_e_v(test_span,i)= abs(d1_real - d1_e);
            end
        end
        RMSE_RML(sample_index,test_span)= sqrt(sum_RML/times);
        RMSE_ML(sample_index,test_span)= sqrt(sum_ML/times);
        CRLB(sample_index,test_span)=sqrt(CRLB_NLoS(R_pd,p,h0,hn,d0_real,d1_real,r_1,m,alpha,beta,test_samples(sample_index)));
    end
    dist_e(sample_index,:)=mean(d1_e_v');
    dist_r(sample_index,:)=mean(d1_r_v');
end

save("figure4.mat","CRLB","RMSE_ML","RMSE_RML","span","dist_r","dist_e");
toc
