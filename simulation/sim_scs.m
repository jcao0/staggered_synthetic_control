
clear
tic 
rng(7)
restoredefaultpath
addpath('functions','data');
alpha_sig = .05;
% warning('off','all')



%% DATA CLEANING
% data = readtable('simulated_data.csv');
data = load('sim_data.mat');
field_names = fieldnames(data);
% data = data(:, {'id','time','D','Y'});

%% PREPARATION
% field_names = field_names{1};
% Initialize result structs to store each NT setting
ATT_struct = struct();
ATT_event_struct = struct();

%% LOOP
for f = 1:numel(field_names)

    field = field_names{f};


    data_field = data.(field);   % each data_field contains one N-T combination
    
    num_sims = size(data_field, 3);     % # of simulation
    fprintf('Processing case: %s, nsims = %d\n', field, num_sims);

    % Prepare storage for this field
%     ATT = cell(num_sims,1);
    ATT_event = cell(num_sims,1);
    
    %%PARALLEL LOOP OVER SIMULATIONS
    parfor i = 1:num_sims
        t0 = tic

        % Print progress every 200 simulations
        if mod(i,200) == 0
            disp(['  Simulation ', num2str(num_sims-i), '/', num2str(num_sims), ' completed...']);
        end

        %%Extract simulation i
        data_i = data_field(:,:,i);
        
        eff = (1:10)';  % true event-tim effect 
        att = unique(data_i(:,5));  % true ATT
        %%RESHAPE DATA
        N = length(unique(data_i(:,1)));
        T1 = length(unique(data_i(:,2)));
        Y = reshape(data_i(:,3),T1,N)'; 
        D = reshape(data_i(:,4),T1,N)'; 
        
        %%DEFINE TREATMENT PERIODS
        T = find(sum(D,1),1)-1;
        S_max = T1-T;
        S = S_max;
        T1 = T + S;
        Y = Y(:,1:T1);
        D = D(:,1:T1);
        D_S = D(:,T+1:T+S);
        
        %%ESTIMATION
%         output = att_ci(Y,D,S,alpha_sig);
%         att_hat = output.att_hat;
%         
        output_t = att_event_ci(Y,D,S,alpha_sig);
        att_hat_t = output_t.att_hat;
        
        %%CONFIDENCE INTERVALS
%         ci_l = att_hat - output.ub;
%         ci_u = att_hat - output.lb;
        ci_l_t = att_hat_t - output_t.ub;
        ci_u_t = att_hat_t - output_t.lb;
        
        %%COLLECT RESULTS
%       ATT{i} = [att, att_hat, att_hat-att, (att_hat-att).^2, ci_l, ci_u, (ci_l<=att) & (ci_u>=att)];
%       ATT_event{i} = [(1:S)', eff(1:S), att_hat_t, att_hat_t - eff(1:S), (att_hat_t-eff(1:S)).^2, ci_l_t, ci_u_t, (ci_l_t<=eff(1:S)) & (ci_u_t>=eff(1:S))];
%       ATT{i} = [att_hat, att_hat-att, (att_hat-att).^2, (ci_l<=att) & (ci_u>=att)];
        ATT_event{i} = [(1:S)', att_hat_t - eff(1:S), att_hat_t, (att_hat_t-eff(1:S)).^2, (ci_l_t<=eff(1:S)) & (ci_u_t>=eff(1:S))];
        
        elapsed_time(i,1) = toc(t0);  
        fprintf('Sim %d | elapsed %.2f sec\n', ...
            i,  elapsed_time(i,1));
    end
    
    %%STORE RESULTS IN STRUCT
%     ATT_struct.(field) = ATT;
    ATT_event_struct.(field) = ATT_event;
end
disp('Finished')
%%
save('sim_result_SSC.mat', 'ATT_struct', 'ATT_event_struct')

%% 
sim_mat = fieldnames(ATT_event_struct);
sim_mat = sim_mat{1};
sim_mat = ATT_event_struct.(sim_mat)(cellfun(@(x) size(x,1) == 10, ATT_event_struct.(sim_mat)));
sim_mat = cat(3, sim_mat{:});

ATT = squeeze(sim_mat(:, 3, :));   % N × sims
ATT_Var = var(ATT, 0, 2); 


elapsed_time = repmat(sum(elapsed_time), length(ATT_Var), 1);
%%
a = mean(sim_mat, 3);
a(:,3) = ATT_Var;    
a = [a, elapsed_time];
a = num2cell(a);
a = cell2table(a, 'VariableNames', ...
    {'event_time','bias','var', 'se','coverage_rate', 'time_total'});
%%
a = ATT_event_struct.N50_T050(cellfun(@(x) size(x,1) == 10, ATT_event_struct.N50_T050));
a = cat(3, a{:});
a = mean(a, 3);

%% test
data = data.N15_T0100(:,:,3);
N = length(unique(data(:,1))); % number of units
T1 = length(unique(data(:,2))); % number of times
Y = reshape(data(:,3),T1,N)'; 
D = reshape(data(:,4),T1,N)'; 
%%
T = find(sum(D),1)-1; % number of pre-treatment periods
S_max = T1-T; % maximum number of post-treatment periods
% S = min(S_max,round(T/2)); % truncate data to avoid extrapolating too far
S = S_max;
T1 = T+S; % total number of periods
Y = Y(:,1:T1); % outcome 
D = D(:,1:T1); % all-time treatment status
D_S = D(:,T+1:T+S); % post-treatment treatment status
%% number of treated units
a = sum(D_S,2);
a = length(find(a));

%% ESTIMATION


output = att_ci(Y,D,S,alpha_sig);
att_hat = output.att_hat;
ub = output.ub;
lb = output.lb;


output_t = att_event_ci(Y,D,S,alpha_sig);
att_hat_t = output_t.att_hat;
ub_t = output_t.ub;
lb_t = output_t.lb;



%% OUTPUT

% % treatment effects estimates (starting from T+1)
% te_mat_hat 

% event time ATT estimates and confidence intervals
disp('[event time, att estimate, confidence interval]')
disp([(1:S)',att_hat_t,att_hat_t-ub,att_hat_t-lb])


% ATT estimates and confidence intervals
disp('[att estimate, confidence interval]')
disp([att_hat,att_hat-ub,att_hat-lb])




%% event time ATT plot
figure
p1 = errorbar(1:S,att_hat_t,ub_t,-lb_t,'CapSize',10);
hold on 
plot(-10:0.1:10,-10:0.1:10,'--k')
plot([.5,S+.5],[0,0],'--k')
hold off
p1.Marker = 'o';
p1.MarkerSize = 6;
p1.LineWidth = 2;
xlabel('event time','FontSize',15)
ylabel('ATT estimates','FontSize',15)
xlim([.5 S+.5])
%% ATT plot
%orginal results in the orginal paper
beta = 0.941;
sd = 0.528;
z = 1.645;   % 90% CI
ci_low = beta - z*sd;
ci_high = beta + z*sd;


figure
hold on 
p1 = errorbar(S/3,att_hat,ub,-lb,'CapSize',10);
p2 = errorbar(2*S/3,beta,z*sd,'CapSize',10);
plot([.5,S+.5],[0,0],'--k')
hold off
p1.Marker = 'o';
p1.MarkerSize = 6;
p1.LineWidth = 2;
p2.Marker = 'o';
p2.MarkerSize = 6;
p2.LineWidth = 2;
ylabel('ATT estimates','FontSize',15)
ax = gca;
ax.XColor = 'none';    % makes the x-axis line invisible
legend([p1 p2], {'SSC', 'Original:General synthetic control (Xu, 2017)'}, ...
       'Location', 'northwest', ...
       'FontSize', 8);


%% residuals/treatment effects plot

% residuals/treatment effects of ever-treated units plot
y_min = min(min(res_mat));
y_max = max(max(res_mat));

figure
subplot(2,1,1)
p = plot(res_mat','LineWidth',2);
hold on 
plot([0,T+S+1],[0,0],'--k')
vline(T,'--k');
hold off
xlim([1,T+S])
ylim([y_min,y_max])
xlabel('time','FontSize',15)
ylabel('treatment effects','FontSize',15)
title('Calendar Time','FontSize',15)


% residuals/treatment effects of ever-treated units plot (recentered at 
% first treatment time)
subplot(2,1,2)
x_min = min(min(event_time_mat));
x_max = max(max(event_time_mat));
plot([x_min-1,x_max+1],[0,0],'--k')
hold on
for i = 1 : size(res_mat,1)
    plot(event_time_mat(i,:),res_mat(i,:),'LineWidth',2)
end
vline(0,'--k');
hold off
x_min_plot = -5; % left limit of plot 
xlim([x_min_plot, x_max])
ylim([y_min,y_max])
xlabel('time relative to first treated','FontSize',15)
ylabel('treament effects','FontSize',15)
title('Event Time','FontSize',15)



%% results by other methods
gsyn = readtable('gsyn.csv');
gsyn.ub = gsyn.ATT - gsyn.CI_lower;
gsyn.lb = gsyn.ATT - gsyn.CI_upper;

asyn = readtable('asyn.csv');
asyn = asyn(asyn.Time>=0,:);
asyn.ub = asyn.Estimate - asyn.lower_bound;
asyn.lb = asyn.Estimate - asyn.upper_bound;

did = readtable('did.csv');
did.ub = did.att - did.CI_lower;
did.lb = did.att - did.CI_upper;

%% plot1 - absolute level treatment effects for two policies
figure(1)
hold on

% --- SSC (coral-red) ---
p1 = errorbar(1:S, att_hat_t, ub_t, -lb_t, ...
    'Color',[1,0.4,0.3],'CapSize',10); % errorbar 1
p1.Marker = 'o'; p1.LineStyle = 'none'; p1.LineWidth = 2;

p2 = plot(1:S, att_hat_t, 'Color',[1,0.4,0.3],'LineWidth',2); % line 2

% --- GSC (teal) ---
p3 = errorbar(1:S, gsyn.ATT, gsyn.ub, -gsyn.lb, ...
    'Color',[0,0.6,0.6],'CapSize',10); % errorbar 3
p3.Marker = 'o'; p3.LineStyle = 'none'; p3.LineWidth = 2;

p4 = plot(1:S, gsyn.ATT, '--', 'Color',[0,0.6,0.6],'LineWidth',2); % line 4

% --- ASC (green) ---
p5 = errorbar(1:size(asyn,1), asyn.Estimate, asyn.ub, -asyn.lb, ...
    'Color',[0.3,0.6,0.2],'CapSize',10); % errorbar 5
p5.Marker = 'o'; p5.LineStyle = 'none'; p5.LineWidth = 2;

p6 = plot(1:size(asyn,1), asyn.Estimate, ':', 'Color',[0.3,0.6,0.2],'LineWidth',2); % line 6

% --- DID (purple) ---
p7 = errorbar(1:S, did.att, did.ub, -did.lb, ...
    'Color',[0.6,0.4,0.8],'CapSize',10); % errorbar 7
p7.Marker = 's'; p7.LineStyle = 'none'; p7.LineWidth = 2;

p8 = plot(1:S, did.att, '-.', 'Color',[0.6,0.4,0.8],'LineWidth',2); % line 8

% --- Axes & labels ---
xlim([0.5, S+0.5])
xlabel('event time','FontSize',15)
ylabel('ATT estimates','FontSize',15)

% --- Legend (reference only line plots) ---
lgd = legend([p2, p4, p6, p8], ...
    {'Staggered synthetic control', ...
     'Staggered diff-in-diff (Callaway and SantAnna, 2018)', ...
     'Augmented synthetic control', ...
     'DID-ASC'}, ...
    'Location','northwest');
lgd.FontSize = 12;

hold off

%%
% save plots
saveas(figure(1),'temp/att_comparison.png')











