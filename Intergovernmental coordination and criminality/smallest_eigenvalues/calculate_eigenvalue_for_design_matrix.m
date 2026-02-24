
% This Matlab script replicates part of results for the empirical replication of:
% "Increasing Intergovernmental Coordination to Fight Crime: Evidence from Mexico"

% The script calculate the smallest eigenvalue of the empirical analog of the design matrix, whose
% invertibility is required by ssc method

% Note: different time periods available for different outcomes
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clear
tic 
rng(7)
restoredefaultpath
addpath('functions','cleaned_data');
alpha_sig = .05;
% warning('off','all')

%%
results_hr = table();
data = readtable('psrm_crime_data.csv');
outcomes = {'hom_all_rate','hom_ym_rate'};
for i =1:length(outcomes)
    outcome = outcomes{i};
    data_i = data(data.time<253, {'idunico','time','Policial',outcome});
    N = length(unique(table2array(unique(data_i(:,1))))); % number of units
    T1 = length(unique(table2array(unique(data_i(:,2))))); % number of times
    Y = reshape(table2array(data_i(:,4)),T1,N)';
    D = reshape(table2array(data_i(:,3)),T1,N)';
    T = find(sum(D),1)-1; % number of pre-treatment periods
    S_max = T1-T; % maximum number of post-treatment periods
    S = S_max; % number of post-treatment periods used
    T1 = T+S; % total number of periods used
    Y = Y(:,1:T1); % outcome
    D = D(:,1:T1); % all-time treatment status
    D_S = D(:,T+1:T+S); % post-treatment treatment status
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
    A = zeros(N,K,S);
    for k = 1 : K
        A(index_mat(k,2),k,index_mat(k,1)) = 1;
    end
    output_t = att_event_ci(Y,D,S,alpha_sig);
    B_hat = output_t.B_hat;
    M_hat = (eye(N)-B_hat)'*(eye(N)-B_hat);
    AMA = 0;
    for s = 1 : S
        AMA = AMA+A(:,:,s)'*M_hat*A(:,:,s);
    end
%     AMA_inv = inv(AMA);
    e = eig(AMA);
    min_eig = min(e);


    result = table(string(outcome), min_eig, ...
        'VariableNames', {'outcome','min_eig'});
    results_hr = [results_hr; result];
    disp(['Outcome completed: ', outcome])
end

%% theft_rate
results_tr = table();
data = readtable('psrm_crime_data.csv');
data.theft_violent_rate = str2double(data.theft_violent_rate);
data.theft_nonviolent_rate = str2double(data.theft_nonviolent_rate);
outcomes = {'theft_violent_rate','theft_nonviolent_rate'};
for i =1:length(outcomes)
    outcome = outcomes{i};
    data_i = data(data.time>=133, {'idunico','time','Policial',outcome});
    N = length(unique(table2array(unique(data_i(:,1))))); % number of units
    T1 = length(unique(table2array(unique(data_i(:,2))))); % number of times
    Y = reshape(table2array(data_i(:,4)),T1,N)';
    D = reshape(table2array(data_i(:,3)),T1,N)';
    T = find(sum(D),1)-1; % number of pre-treatment periods
    S_max = T1-T; % maximum number of post-treatment periods
    S = S_max;
    T1 = T+S; % total number of periods
    Y = Y(:,1:T1); % outcome
    D = D(:,1:T1); % all-time treatment status
    D_S = D(:,T+1:T+S); % post-treatment treatment status
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
    A = zeros(N,K,S);
    for k = 1 : K
        A(index_mat(k,2),k,index_mat(k,1)) = 1;
    end

    output_t = att_event_ci(Y,D,S,alpha_sig);
    B_hat = output_t.B_hat;
    M_hat = (eye(N)-B_hat)'*(eye(N)-B_hat);
    AMA = 0;
    for s = 1 : S
        AMA = AMA+A(:,:,s)'*M_hat*A(:,:,s);
    end
%     AMA_inv = inv(AMA);
    e = eig(AMA);
    min_eig = min(e);


    result = table(string(outcome), min_eig, ...
        'VariableNames', {'outcome','min_eig'});
    results_tr = [results_tr; result];
    disp(['Outcome completed: ', outcome])
end

%% cartel presence
results_cp = table();
data = readtable('psrm_cartel_data.csv');
outcomes = {'presence_strength','co_num','war'};
for i =1:length(outcomes)
    outcome = outcomes{i};
    data_i = data(:, {'idunico','Year','policial',outcome});
    N = length(unique(table2array(unique(data_i(:,1))))); % number of units
    T1 = length(unique(table2array(unique(data_i(:,2))))); % number of times
    Y = reshape(table2array(data_i(:,4)),T1,N)';
    D = reshape(table2array(data_i(:,3)),T1,N)';
    T = find(sum(D),1)-1; % number of pre-treatment periods
    S_max = T1-T; % maximum number of post-treatment periods
    % S = min(S_max,round(T/2)); % truncate data to avoid extrapolating too far
    S = S_max;
    T1 = T+S; % total number of periods
    Y = Y(:,1:T1); % outcome
    D = D(:,1:T1); % all-time treatment status
    D_S = D(:,T+1:T+S); % post-treatment treatment status
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
    A = zeros(N,K,S);
    for k = 1 : K
        A(index_mat(k,2),k,index_mat(k,1)) = 1;
    end
    output_t = att_event_ci(Y,D,S,alpha_sig);
    B_hat = output_t.B_hat;
    M_hat = (eye(N)-B_hat)'*(eye(N)-B_hat);
    AMA = 0;
    for s = 1 : S
        AMA = AMA+A(:,:,s)'*M_hat*A(:,:,s);
    end
%     AMA_inv = inv(AMA);
    e = eig(AMA);
    min_eig = min(e);


    result = table(string(outcome), min_eig, ...
        'VariableNames', {'outcome','min_eig'});
    results_cp = [results_cp; result];
    disp(['Outcome completed: ', outcome])
end
%% collect and save results
results_ssc= [results_hr; results_tr; results_cp];
if ~exist('output', 'dir')
    mkdir('output')
end
writetable(results_ssc,'output/Table1_smallest_eigenvalue.csv')
fprintf('Results saved in /output\n');