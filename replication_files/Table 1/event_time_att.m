tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
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
writetable(weight_table,'table1.csv')








