function d_n = ML_NLoS(ro, K_n, sigma_T_2, q, B, I_1, R, p, A, G, m, q_z, wn_z, r_n, mu_n,d0)
% PARAMETERS
% ro =
% K_n =
% sigma_T_2 =
% q =
% B =
% I_1 =
% R =
% p =
% A =
% G =
% m =
% q_z =
% wn_z =
% r_n =
% mu_n =
% RETURNS
% d_n =

alpha= (sigma_T_2) + (2*q*I_1*B);
beta = (2*q*B);

ups = (R*p*A*G*(m+1)*(q_z^(m+1)))/((d0^(m+3))*(2*pi));
kappa = (R*p*ro*A*G*(m+1)*((q_z - wn_z)^m)*wn_z)/(2*pi*(r_n^m));

somma=sum(mu_n);
somma2=sum(mu_n.*mu_n);


%syms d
%eqn=(somma2 - (2*somma*kappa)/(d*((d + r_n)^2)) + (K_n*kappa^2)/((d^2)*(d + r_n)^4))/(alpha + beta*(kappa/(d*((d + r_n)^2)) + ups)) + log10(alpha + beta*(kappa/(d*((d + r_n)^2)) + ups));
fun = @(d)ML_loglikelihood(d,somma2,somma,kappa,r_n,K_n,alpha,beta,ups);

%1o step: stima distanza utilizzando una prima griglia coarse
Np = 10000; %numero di punti nella griglia
dist_max = 15; %massima distanza per la griglia
span_dist = linspace(0,dist_max,Np); %creo i vari Np punti da testare in maniera equidistanziati tra [0,d_max]
     
ML_dist = zeros(length(span_dist),1); %inizializzo vettore per salvare i valori della log-likelihood per ogni d di tentativo nella griglia



for j=1:length(span_dist)
    
  ML_dist(j) =  fun(span_dist(j)); % eqn(span_dist(j));

end

minimum_ML_dist= min(ML_dist); %trovo il minimo della log likelihood
       
ind_ML_dist = find(ML_dist==minimum_ML_dist,1,'first'); %trovo indice corrispondente al minimo nel vettore dei valori della log likelihood

estimated_dist_MLcoarse = span_dist(ind_ML_dist); %recupero il valore di d che minimizzava, quindi l'output dello stimatore

  
%2o step: utilizzo estimated_dist_MLcoarse come valore iniziale per un
%approccio iterativo tipo "Gradient-based"
x0 = estimated_dist_MLcoarse; %valore iniziale, partiamo dalla guess della stima coarse basata su griglia

%rescale = ML_loglikelihood(x0,...); %quando usi algoritmi iterativi, è spesso utile riscalare, questa poi te la spiego a voce
%fun = @(x)ML_loglikelihood(x,...); %uso un richiamo parametrico in x (si chiamano funzioni implicite) alla funzione che ho utilizzato a riga 10 e che mi calcola log likelihood

%options = optimset('TolFun', 1e-3,'TolX',1e-3,'MaxFunEvals',3000,'display','off');
options = optimoptions('fmincon','FunctionTolerance',1e-12,'StepTolerance',1e-12,'OptimalityTolerance',1e-12);
dist_refined = fmincon(fun,x0,[],[],[],[],0,dist_max,[],options); %stima finale della distanza raffinata tramite ottimizzazione iterativa

d_n=dist_refined;