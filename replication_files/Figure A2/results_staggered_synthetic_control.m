
clear
tic 
rng(7)
restoredefaultpath
addpath('functions');
alpha_sig = .05;
% warning('off','all')



%% STAGGERED SYNTHETIC CONTROL
% DATA CLEANING

data = readtable('data_boardgendereige.csv');
head(data,5)
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


%% READ CS DID RESULTS
cs = readtable('results_staggered_did.csv');
cs = cs(ismember(cs.event_time, [0:S-1]), :);
cs.ub = cs.ci_upper - cs.att;
cs.lb = cs.att - cs.ci_lower;

%% plot - treatment effects from two methods
figure(1)

hold on 
% plot([.5,S+.5],[0,0],'--k');

p1 = errorbar(1:S,att_hat,ub,-lb,...
    'Color',[1,0.4,0.3],'CapSize',10);
p1.Marker = 'o';
p1.LineStyle = 'none';
p1.LineWidth = 2;
p2 = errorbar(1:S,cs.att,cs.lb,cs.ub,...
    'Color',[0,.6,.6],'CapSize',10);
p2.Marker = 'o';
p2.LineStyle = 'none';
p2.LineWidth = 2;
p3 = plot(1:S,att_hat,'Color',[1,0.4,0.3],'LineWidth',2);
p4 = plot(1:S,cs.att,'--','Color',[0,.6,.6],'LineWidth',2);

xlim([.5 S+.5])
xlabel('event time','FontSize',15)
ylabel('ATT estimates','FontSize',15)


% title('Comparison between Staggered Synthetic Control and Staggered Diff-in-diff (Callaway, B. and Sant'Anna, P. H. C., 2021)','FontSize',15)
lgd = legend([p3(1) p4(1)],"Staggered synthetic control","Staggered diff-in-diff (Callaway and Sant’Anna, 2021)",...
    'Location','northwest');
lgd.FontSize = 12;

hold off

% save plots
saveas(figure(1),'FigureA2.png')











