function [a_hat,B_hat] = synthetic_control_batch(Y)
% SYNTHETIC_CONTROL_BATCH calculate all synthetic control weights, using
% each row as the treated and the others as the controls, separately. 

[N,T] = size(Y);
a_hat = zeros(N,1);
B_hat = zeros(N);

for i = 1 : N
    Y_treated = Y(i,:)';
    temp = Y;
    temp(i,:) = [];
    Y_untreated = temp';

    Y_demeaned = Y_treated-mean(Y_treated);
    X_demeaned = Y_untreated-repmat(mean(Y_untreated),T,1);

    b_initial = ones(N-1,1)/(N-1);
    Q = @(b)sum((Y_demeaned-X_demeaned*b).^2); % criterion

    % constraints
    A_eq = ones(1,N-1);
    B_eq = 1;
    LB = zeros(N-1,1);

    options = optimoptions('fmincon','Display','none');
    b_hat = fmincon(Q,b_initial,[],[],A_eq,B_eq,LB,[],[],options);

    a_hat(i) = mean(Y_treated)-mean(Y_untreated)*b_hat;
    b_hat = [b_hat(1:i-1);0;b_hat(i:end)];
    B_hat(i,:) = b_hat';
end