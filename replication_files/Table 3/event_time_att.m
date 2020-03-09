
clear
tic 
rng(7)
restoredefaultpath
addpath('functions');
alpha_sig = .05;
% warning('off','all')


%% DATA CLEANING

data = readtable('data_boardgendereige.csv');
% head(data,5)
N = length(unique(table2array(unique(data(:,1))))); % number of units
T1 = length(unique(table2array(unique(data(:,3))))); % number of units
Y = reshape(table2array(data(:,7)),T1,N)'; 
D = reshape(table2array(data(:,4)),T1,N)'; 
T = find(sum(D),1)-1; % number of pre-treatment periods
S_max = T1-T; % maximum number of post-treatment periods

S = S_max-3; % truncate data to avoid extrapolating too far
T1 = T+S; % total number of periods
Y = Y(:,1:T1); % outcome 
D = D(:,1:T1); % all-time treatment status
D_S = D(:,T+1:T+S); % post-treatment treatment status


%% ESTIMATION

output = att_event_ci(Y,D,S,alpha_sig);
te_mat_hat = output.te_mat_hat;
att_hat = output.att_hat;
ub = output.ub;
lb = output.lb;
res_mat = output.res_mat;
event_time_mat = output.event_time_mat;
B_hat = output.B_hat;

% display weights for a certian unit
id_unit = 9; % unit of interest
weights = round(B_hat(id_unit,:)',4);
units = unique(table2array(unique(data(:,2))));
weight_table = table(units,weights);
weight_table(id_unit,:) = []; % delete itself
weight_table
writetable(weight_table,'table3.csv')

%% OUTPUT

% % treatment effects estimates (starting from T+1)
% te_mat_hat 

% event time ATT estimates and confidence intervals
disp('[event time, att estimate, confidence interval]')
disp([(1:S)',att_hat,att_hat-ub,att_hat-lb])

%% event time ATT plot
figure
p1 = errorbar(1:S,att_hat,ub,-lb);
hold on 
plot([.5,S+.5],[0,0],'--k')
hold off
p1.Marker = 'o';
p1.MarkerSize = 6;
p1.LineWidth = 1.5;
xlabel('event time')
ylabel('ATT estimates')
xlim([.5 S+.5])

%% residuals/treatment effects plot

% residuals/treatment effects of ever-treated units plot
y_min = min(min(res_mat));
y_max = max(max(res_mat));

figure
subplot(2,1,1)
p = plot(res_mat');
hold on 
plot([0,T+S+1],[0,0],'--k')
vline(T,'--k');
hold off
xlim([1,T+S])
ylim([y_min,y_max])
xlabel('time')
ylabel('treatment effects')
title('Calendar Time')


% residuals/treatment effects of ever-treated units plot (recentered at 
% first treatment time)
subplot(2,1,2)
x_min = min(min(event_time_mat));
x_max = max(max(event_time_mat));
plot([x_min-1,x_max+1],[0,0],'--k')
hold on
for i = 1 : size(res_mat,1)
    plot(event_time_mat(i,:),res_mat(i,:))
end
vline(0,'--k');
hold off
x_min_plot = -5; % left limit of plot 
xlim([x_min_plot, x_max])
ylim([y_min,y_max])
xlabel('time relative to first treated')
ylabel('treament effects')
title('Event Time')




toc









