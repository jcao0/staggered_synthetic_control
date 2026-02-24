% This Matlab script replicates part of results for the empirical replication of:
% "Increasing Intergovernmental Coordination to Fight Crime: Evidence from Mexico"

% The script includes results estimated by Staggered Synthetic Control

% Note: different time periods available for different outcomes
clear
tic 
rng(7)
restoredefaultpath
addpath('functions','data/processed/');
alpha_sig = .05;
% warning('off','all')


%% homicide_rate
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
    output_t = att_event_ci(Y,D,S,alpha_sig);
    te_mat_hat = output_t.te_mat_hat;
    att_hat_t = output_t.att_hat;
    ub_t = output_t.ub;
    lb_t = output_t.lb;
    outcome_col = repmat(string(outcome), S, 1);
    TN_col      = repmat(T/N, S, 1);
    N_col       = repmat(N, S, 1);
    T_col       = repmat(T, S, 1);
    S_col       = repmat(S, S, 1);
    event_time  = (1:S)';
    result = table(outcome_col, TN_col, N_col, T_col, S_col, event_time, att_hat_t, att_hat_t-ub_t, att_hat_t-lb_t, ...
          'VariableNames', {'outcome','T/N','N','T','S','event time','att estimate', 'confidence interval_l', 'confidence interval_u'});
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
    % S = min(S_max,round(T/2)); % truncate data to avoid extrapolating too far
    S = S_max;
    T1 = T+S; % total number of periods
    Y = Y(:,1:T1); % outcome
    D = D(:,1:T1); % all-time treatment status
    D_S = D(:,T+1:T+S); % post-treatment treatment status
    output_t = att_event_ci(Y,D,S,alpha_sig);
    te_mat_hat = output_t.te_mat_hat;
    att_hat_t = output_t.att_hat;
    ub_t = output_t.ub;
    lb_t = output_t.lb;
    outcome_col = repmat(string(outcome), S, 1);
    TN_col      = repmat(T/N, S, 1);
    N_col       = repmat(N, S, 1);
    T_col       = repmat(T, S, 1);
    S_col       = repmat(S, S, 1);
    event_time  = (1:S)';
    result = table(outcome_col, TN_col, N_col, T_col, S_col, event_time, att_hat_t, att_hat_t-ub_t, att_hat_t-lb_t, ...
          'VariableNames', {'outcome','T/N','N','T','S','event time','att estimate', 'confidence interval_l', 'confidence interval_u'});
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
    output_t = att_event_ci(Y,D,S,alpha_sig);
    te_mat_hat = output_t.te_mat_hat;
    att_hat_t = output_t.att_hat;
    ub_t = output_t.ub;
    lb_t = output_t.lb;
    outcome_col = repmat(string(outcome), S, 1);
    TN_col      = repmat(T/N, S, 1);
    N_col       = repmat(N, S, 1);
    T_col       = repmat(T, S, 1);
    S_col       = repmat(S, S, 1);
    event_time  = (1:S)';
    result = table(outcome_col, TN_col, N_col, T_col, S_col, event_time, att_hat_t, att_hat_t-ub_t, att_hat_t-lb_t, ...
          'VariableNames', {'outcome','T/N','N','T','S','event time','att estimate', 'confidence interval_l', 'confidence interval_u'});
    results_cp = [results_cp; result];
    disp(['Outcome completed: ', outcome])
end

%% collect and save results
results_ssc= [results_hr; results_tr; results_cp];

writetable(results_ssc,'results/results_ssc.csv')
