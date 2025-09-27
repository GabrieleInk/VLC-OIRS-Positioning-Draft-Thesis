%% Performance di localizzazione basata sulla legenda di colori (freddo (~1mm) - caldo (~4mm) quindi errore commettibile più alto) con PEB, flusso luminoso fisso a 1000lm  e 2 iterazioni dell'algoritmo di adaptive beam steering


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

granularity=100; % granularity=1000;

res=x_max/granularity;
x_scan=0:res:x_max;
y_scan=0:res:y_max;
[x,y]=meshgrid(0:res:x_max,0:res:y_max);

z_ris=1.35;
theta1 = -30;  % rotazione degli specchietti di sinistra sul piano (se presenti)
theta2 = 30;  % rotazione degli specchietti di destra sul piano (se presenti)

for RIS_config=1:3 % for RIS_config=1:1 per eseguire solo il caso con 4 OIRS
    clear RIS Norm alpha beta w h;

    % 4 OIRS
    if RIS_config==1
        RIS=zeros(4,3);
        Norm=zeros(4,3);
        alpha=zeros(4,1);
        beta=zeros(4,1);
        w=zeros(4,1);
        h=zeros(4,1);

        RIS(1,:)    = [x_max/2      , 0            , z_ris];
        Norm(1,:)=[0 1 0];
        alpha(1) = 0;
        beta(1) = 0;
        w(1)=1;
        h(1)=1;

        RIS(2,:)    = [0            , y_max/2      , z_ris];
        Norm(2,:)=[1 0 0];
        alpha(2) = 0;
        beta(2) = 0;
        w(2)=1;
        h(2)=1;

        RIS(3,:)    = [x_max        , y_max/2      , z_ris];
        Norm(3,:)=[-1 0 0];
        alpha(3) = 0;
        beta(3) = 0;
        w(3)=1;
        h(3)=1;

        RIS(4,:)    = [x_max/2      , y_max        , z_ris];
        Norm(4,:)=[0 -1 0];
        alpha(4) = 0;
        beta(4) = 0;
        w(4)=1;
        h(4)=1;
    
    % 8 OIRS    
    elseif RIS_config==2
        RIS=zeros(8,3);
        Norm=zeros(8,3);
        alpha=zeros(8,1);
        beta=zeros(8,1);
        w=zeros(8,1);
        h=zeros(8,1);

        RIS(1,:)    = [(x_max/5)*3      , 0            , z_ris];
        Norm(1,:)=[0 1 0];
        alpha(1) = 0;
        beta(1) = 0;
        w(1)=1;
        h(1)=1;

        RIS(5,:)    = [(x_max/5)*2      , 0            , z_ris];
        Norm(5,:)=[0 1 0];
        alpha(5) = -0;
        beta(5) = 0;
        w(5)=1;
        h(5)=1;

        RIS(2,:)    = [0            , (y_max/5)*3      , z_ris];
        Norm(2,:)=[1 0 0];
        alpha(2) = -0;
        beta(2) = 0;
        w(2)=1;
        h(2)=1;

        RIS(6,:)    = [0            , (y_max/5)*2      , z_ris];
        Norm(6,:)=[1 0 0];
        alpha(6) = 0;
        beta(6) = 0;
        w(6)=1;
        h(6)=1;

        RIS(3,:)    = [x_max        , (y_max/5)*3      , z_ris];
        Norm(3,:)=[-1 0 0];
        alpha(3) = 0;
        beta(3) = 0;
        w(3)=1;
        h(3)=1;

        RIS(7,:)    = [x_max        , (y_max/5)*2      , z_ris];
        Norm(7,:)=[-1 0 0];
        alpha(7) = -0;
        beta(7) = 0;
        w(7)=1;
        h(7)=1;

        RIS(4,:)    = [(x_max/5)*3      , y_max        , z_ris];
        Norm(4,:)=[0 -1 0];
        alpha(4) = -0;
        beta(4) = 0;
        w(4)=1;
        h(4)=1;

        RIS(8,:)    = [(x_max/5)*2      , y_max        , z_ris];
        Norm(8,:)=[0 -1 0];
        alpha(8) = 0;
        beta(8) = 0;
        w(8)=1;
        h(8)=1;

    % 12 OIRS    
    elseif RIS_config==3
        RIS=zeros(12,3);
        Norm=zeros(12,3);
        alpha=zeros(12,1);
        beta=zeros(12,1);
        w=zeros(12,1);
        h=zeros(12,1);

        RIS(1,:)    = [(x_max/10)*7      , 0            , z_ris];
        Norm(1,:)=[0 1 0];
        alpha(1) = 0;
        beta(1) = 0;
        w(1)=1;
        h(1)=1;

        RIS(5,:)    = [(x_max/10)*5      , 0            , z_ris];
        Norm(5,:)=[0 1 0];
        alpha(5) = -0;
        beta(5) = 0;
        w(5)=1;
        h(5)=1;

        RIS(9,:)    = [(x_max/10)*3      , 0            , z_ris];
        Norm(9,:)=[0 1 0];
        alpha(9) = -0;
        beta(9) = 0;
        w(9)=1;
        h(9)=1;

        RIS(2,:)    = [0            , (y_max/10)*7      , z_ris];
        Norm(2,:)=[1 0 0];
        alpha(2) = -0;
        beta(2) = 0;
        w(2)=1;
        h(2)=1;

        RIS(6,:)    = [0            , (y_max/10)*5      , z_ris];
        Norm(6,:)=[1 0 0];
        alpha(6) = 0;
        beta(6) = 0;
        w(6)=1;
        h(6)=1;

        RIS(10,:)    = [0            , (y_max/10)*3      , z_ris];
        Norm(10,:)=[1 0 0];
        alpha(10) = 0;
        beta(10) = 0;
        w(10)=1;
        h(10)=1;

        RIS(3,:)    = [x_max        , (y_max/10)*7      , z_ris];
        Norm(3,:)=[-1 0 0];
        alpha(3) = 0;
        beta(3) = 0;
        w(3)=1;
        h(3)=1;

        RIS(7,:)    = [x_max        , (y_max/10)*5      , z_ris];
        Norm(7,:)=[-1 0 0];
        alpha(7) = -0;
        beta(7) = 0;
        w(7)=1;
        h(7)=1;

        RIS(11,:)    = [x_max        , (y_max/10)*3      , z_ris];
        Norm(11,:)=[-1 0 0];
        alpha(11) = -0;
        beta(11) = 0;
        w(11)=1;
        h(11)=1;

        RIS(4,:)    = [(x_max/10)*7      , y_max        , z_ris];
        Norm(4,:)=[0 -1 0];
        alpha(4) = -0;
        beta(4) = 0;
        w(4)=1;
        h(4)=1;

        RIS(8,:)    = [(x_max/10)*5      , y_max        , z_ris];
        Norm(8,:)=[0 -1 0];
        alpha(8) = 0;
        beta(8) = 0;
        w(8)=1;
        h(8)=1;

        RIS(12,:)    = [(x_max/10)*3      , y_max        , z_ris];
        Norm(12,:)=[0 -1 0];
        alpha(12) = 0;
        beta(12) = 0;
        w(12)=1;
        h(12)=1;
    end
  
    [alp, bet] = noiseParameters(q_0, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);

    result(RIS_config).crlb = generatePEB(x_scan, y_scan, LED, RIS, Norm, Phi_FoV, a, rho, Psi, A_pd, T_of,R_pd,p,alp,bet,K0,KN,m,alpha,beta,w,h);
                              
end

save("figure7_1.mat","result","x","y");

elapsedtime7_1=toc
if isfile("savings.mat")
    load("savings.mat");  % carica la struct "savings"
else
    savings = struct();   % crea la struct se non esiste
end

savings.Figure7_1_generate = elapsedtime7_1;  % aggiungi/aggiorna campo

save("savings.mat","savings");