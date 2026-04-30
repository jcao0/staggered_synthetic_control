function [output] = ssc(Y,D)

[N,T1] = size(D);
T = find(sum(D),1)-1;
S = T1-T;

Y_T = Y(:,1:T);
Y_S = Y(:,T+1:end);
D_S = D(:,T+1:end);

% index matrix - [time,unit,event_time]
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
tau_hat = temp1\temp2; % heterogeneous treatment effects

u_hat = Y_T-(repmat(a_hat,1,T)+B_hat*Y_T); % residuals

V_mat = zeros(K,T-S);
for t = 1 : T-S
    temp2 = 0;
    for s = 1 : S
        temp2 = temp2+A(:,:,s)'*(eye(N)-B_hat)'*u_hat(:,t+s);
    end
    V_mat(:,t) = temp1\temp2;
end

output.tau_hat = tau_hat;
output.V_mat = V_mat;
output.B_hat = B_hat;










