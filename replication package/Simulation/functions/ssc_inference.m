function [gamma_hat,lb,ub] = ssc_inference(L,alpha,tau_hat,V_mat)

gamma_hat = L*tau_hat;
null_distr = L*V_mat;

lb = quantile(null_distr,alpha/2);
ub = quantile(null_distr,1-alpha/2);















