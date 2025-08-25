function loglike=ML_loglikelihood(d,somma2,somma,kappa,r_n,K_n,alpha,beta,ups)


loglike=(somma2 - (2*somma*kappa)/(d*((d + r_n)^2)) + (K_n*kappa^2)/((d^2)*(d + r_n)^4))/(alpha + beta*(kappa/(d*((d + r_n)^2)) + ups)) + log10(alpha + beta*(kappa/(d*((d + r_n)^2)) + ups));