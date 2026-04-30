function [output] = att_ci(Y,D,S,alpha_sig)

[N,~] = size(D);
T = find(sum(D),1)-1;

% truncate data at T+S
Y = Y(:,1:T+S);
D = D(:,1:T+S);

Y_T = Y(:,1:T);
Y_S = Y(:,T+1:end);
D_S = D(:,T+1:end);

% index matrix - [post_time,unit,event_time]
K = sum(sum(D)); % dimension of gamma
index_mat = zeros(K,3);
ind = 0;
for s = 1 : S
    for i = 1 : N
        if D_S(i,s)==1
            ind = ind+1;
            index_mat(ind,:) = [s,i,sum(D_S(i,1:s))];
        end
    end
end

% treatment structure - A_s matrices 
A = zeros(N,K,S);
for k = 1 : K
    A(index_mat(k,2),k,index_mat(k,1)) = 1;
end

% synthetic control weights
[a_hat,B_hat] = synthetic_control_batch(Y_T);
M_hat = (eye(N)-B_hat)'*(eye(N)-B_hat);

% estimation
temp1 = 0;
for s = 1 : S
    temp1 = temp1+A(:,:,s)'*M_hat*A(:,:,s);
end

temp2 = 0;
for s = 1 : S
    temp2 = temp2+A(:,:,s)'*(eye(N)-B_hat)'*((eye(N)-B_hat)*...
        Y_S(:,s)-a_hat);
end

gamma_hat = temp1\temp2;

% linear transformation from gamma to alpha
L = zeros(1,K);
L(1,:) = 1/K;

att_hat = L*gamma_hat;

u_hat = Y_T-(repmat(a_hat,1,T)+B_hat*Y_T);

V_mat = zeros(K,T-S);

for t = 1 : T-S
    temp2 = 0;
    for s = 1 : S
        temp2 = temp2+A(:,:,s)'*(eye(N)-B_hat)'*u_hat(:,t+s);
    end
    V_mat(:,t) = temp1\temp2;
end

null_distr = L*V_mat;

lb = quantile(null_distr,alpha_sig/2);
ub = quantile(null_distr,1-alpha_sig/2);

null_distr_abs = abs(null_distr);
att_hat_abs = abs(att_hat);
p = sum(null_distr_abs > att_hat_abs) / numel(null_distr_abs);



te_mat_hat = zeros(N,S);
for s = 1 : S
    te_mat_hat(:,s) = A(:,:,s)*gamma_hat;
end

res_mat = [u_hat,te_mat_hat];
res_mat(sum(D,2)==0,:) = [];


output.te_mat_hat = te_mat_hat;
output.att_hat = att_hat;
output.ub = ub;
output.lb = lb;
output.res_mat = res_mat;
output.B_hat = B_hat;
output.p = p;











