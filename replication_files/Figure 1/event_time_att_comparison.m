tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clear
tic 
rng(7)
addpath(genpath('functions'));
alpha_sig = .05;
% warning('off','all')


%% DATA CLEARNING

data = readtable('data_boardgendereige.csv');
% head(data,5)
N = length(unique(table2array(unique(data(:,1))))); % number of units
T1 = length(unique(table2array(unique(data(:,3))))); % number of periods
Y = reshape(table2array(data(:,7)),T1,N)'; 
D = reshape(table2array(data(:,4)),T1,N)'; 
D1 = reshape(table2array(data(:,5)),T1,N)'; 
D2 = reshape(table2array(data(:,6)),T1,N)'; 
T = find(sum(D),1)-1; % number of pre-treatment periods
S_max = T1-T; % maximum number of post-treatment periods

S = S_max-3; % truncate data to avoid extrapolating too far
T1 = T+S; % total number of periods
Y = Y(:,1:T1); % outcome 
D = D(:,1:T1); % all-time treatment status
D_S = D(:,T+1:T+S); % post-treatment treatment status
D1 = D1(:,1:T1);
D2 = D2(:,1:T1);
D1_S = D1(:,T+1:T+S);
D2_S = D2(:,T+1:T+S);

%% ESTIMATION AND INFERENCE

% apply the synthetic control method
output = ssc(Y,D);
tau_hat = output.tau_hat;
V_mat = output.V_mat;
B_hat = output.B_hat;

% display weights for a certian unit
id_unit = 9; % unit of interest
weights = round(B_hat(id_unit,:)',4);
units = unique(table2array(unique(data(:,2))));
weight_table = table(units,weights);
weight_table(id_unit,:) = []; % delete itself
weight_table

% index matrix - [time,unit,event_time,policy]
K = sum(sum(D)); % dimension of tau
index_mat = zeros(K,4);
ind = 0;
for s = 1 : S
    for i = 1 : N
        if D_S(i,s)==1
            ind = ind+1;
            if D1_S(i,s) == 1
                index_mat(ind,:) = [s,i,sum(D_S(i,1:s)),1];
            else
                index_mat(ind,:) = [s,i,sum(D_S(i,1:s)),2];
            end
                
        end
    end
end

% Policy 1
S1 = max(sum(D1_S,2));
gamma_hat_1 = zeros(S1,1);
lb_1 = zeros(S1,1);
ub_1 = zeros(S1,1);
for s = 1 : S1
    ind = (index_mat(:,3) == s).*(index_mat(:,4)==1);
    L = ind'/sum(ind);
    [gamma_hat,lb,ub] = ssc_inference(L,alpha_sig,tau_hat,V_mat);
    gamma_hat_1(s) = gamma_hat;
    lb_1(s) = lb;
    ub_1(s) = ub;
end

% Policy 2
S2 = max(sum(D2_S,2));
gamma_hat_2 = zeros(S2,1);
lb_2 = zeros(S2,1);
ub_2 = zeros(S2,1);
for s = 1 : S2
    ind = (index_mat(:,3) == s).*(index_mat(:,4)==2);
    L = ind'/sum(ind);
    [gamma_hat,lb,ub] = ssc_inference(L,alpha_sig,tau_hat,V_mat);
    gamma_hat_2(s) = gamma_hat;
    lb_2(s) = lb;
    ub_2(s) = ub;
end

% difference
S_min = min(S1,S2);
gamma_hat_diff = zeros(S_min,1);
lb_vec = zeros(S_min,1);
ub_vec = zeros(S_min,1);
for s = 1 : S_min
    ind1 = (index_mat(:,3) == s).*(index_mat(:,4)==1);
    L1 = ind1'/sum(ind1);
    ind2 = (index_mat(:,3) == s).*(index_mat(:,4)==2);
    L2 = ind2'/sum(ind2);
    L = L1-L2;
    [gamma_hat,lb,ub] = ssc_inference(L,alpha_sig,tau_hat,V_mat);
    gamma_hat_diff(s) = gamma_hat;
    lb_vec(s) = lb;
    ub_vec(s) = ub;
end
%% data for plot
output = att_event_ci(Y,D,S,alpha_sig);
te_mat_hat = output.te_mat_hat;
att_hat = output.att_hat;
ub = output.ub;
lb = output.lb;
res_mat = output.res_mat;
event_time_mat = output.event_time_mat;
B_hat = output.B_hat;

D1_temp = D1_S;
D2_temp = D2_S;
ind_temp = sum(D,2)==0;
D1_temp(ind_temp,:) = [];
D2_temp(ind_temp,:) = [];
res_mat_1 = res_mat((sum(D1_temp,2)>0),:);
res_mat_2 = res_mat((sum(D2_temp,2)>0),:);% does not include germany, as it takes both quota and disclosure


% residuals/treatment effects of ever-treated units plot
y_min = min(min(res_mat));
y_max = max(max(res_mat));

% select units of policy 1
units_p1 = (sum(D1_S,2) > 0);
units_p1 = units(units_p1);
% select units of policy 2
units_p2 = (sum(D2_S,2) > 0);
units_p2 = units(units_p2);

% Full country names
countries = {...
    'Austria', 'Belgium', 'Bulgaria', 'Croatia', 'Cyprus', ...
    'Czech Republic', 'Denmark', 'Estonia', 'France', 'Germany', ...
    'Greece', 'Hungary', 'Ireland', 'Italy', 'Latvia', 'Lithuania', ...
    'Luxembourg', 'Malta', 'Netherlands', 'Poland', 'Portugal', ...
    'Romania', 'Slovakia', 'Slovenia', 'Sweden', 'United Kingdom'};

% Corresponding abbreviations (2-letter or ISO codes)
abbr = {...
    'AT', 'BE', 'BG', 'HR', 'CY', ...
    'CZ', 'DK', 'EE', 'FR', 'DE', ...
    'GR', 'HU', 'IE', 'IT', 'LV', 'LT', ...
    'LU', 'MT', 'NL', 'PL', 'PT', ...
    'RO', 'SK', 'SI', 'SE', 'UK'};

% Create mapping
countryMap = containers.Map(countries, abbr);
units_p1 = cellfun(@(x) countryMap(x), units_p1, 'UniformOutput', false);
units_p2 = cellfun(@(x) countryMap(x), units_p2, 'UniformOutput', false);
%% plot  - residual plo
% residuals/treatment effects of ever-treated units plot
fig = figure('Position', [616, 598, 33*70/3, 33*70/4]);
policy = {'Quota', 'Disclosure'};
colors     = {[1, 0.4, 0.3], [0,.6,.6]};
linestyles = {'-', '--'};
markers = {'o','+','^','s'};

subplot(2,1,1)
hold on
plot([0,T+S+1],[0,0],'--k','HandleVisibility','off');

n1 = size(res_mat_1, 1);
n2 = size(res_mat_2, 1);
h = nan(max(n1,n2),length(policy));

% plot for each country
for i = 1:n1
      % Plot
      h(i,1) = plot(res_mat_1(i, :), ...
        '-',...
        'Color', [1, 0.4, 0.3], ...
        'LineWidth', 2, ...
        'Marker', markers{i}, ...
        'MarkerSize', 10);
end
for i = 1:n2
      % Plot
      h(i,2) = plot(res_mat_2(i, :), ...
        '--',...
        'Color', [0,.6,.6], ...
        'LineWidth', 2, ...
        'Marker', markers{i}, ...
        'MarkerSize', 10);
end

xlim([1-0.5,T+S+0.5])
ylim([y_min-0.5,y_max+1])
vline(T,'--k');
xlabel('time','FontSize',15)
ylabel('treatment effects','FontSize',15)
title('Calendar Time','FontSize',15)

h1 = ...
legendflex(h(1:n1,1), units_p1, ...
'ref', gca, ...
'anchor', {'nw','nw'}, ...
'buffer', [5 -5], ...
'fontsize',12, ...
'xscale',1.2, ...
'title', 'Quota');
sz = get(h1, 'position');
wid_box = 1.2*sz(3);
sz(3) = wid_box;     % Make the width larger
set(h1,'position', sz);


h2 = ...
legendflex(h(1:n2,2), units_p2, ...
'ref', h1, ...
'anchor', {'sw','nw'}, ...
'buffer', [0 -8], ...
'fontsize', 12, ...
'xscale',1.2, ...
'box','on',...
'title', 'Disclosure');
sz = get(h2, 'position');
sz(3) = wid_box;     % Make the heigth larger
set(h2,'position', sz);

box on
hold off



% residuals/treatment effects of ever-treated units plot (recentered at
% first treatment time)
event_time_mat1 = event_time_mat((sum(D1_temp,2)>0),:);
event_time_mat2 = event_time_mat((sum(D2_temp,2)>0),:);
res_mat_1 = res_mat((sum(D1_temp,2)>0),:);
res_mat_2 = res_mat((sum(D2_temp,2)>0),:);
x_min = min(min(event_time_mat));
x_max = max(max(event_time_mat));

subplot(2,1,2)
hold on
plot([x_min-1,x_max+1],[0,0],'--k')

n1 = size(event_time_mat1, 1);
n2 = size(event_time_mat2, 1);
h = nan(max(n1,n2),length(policy));


% plot for each country
for i = 1:n1
      % Plot
      h(i,1) = plot(event_time_mat1(i, :),res_mat_1(i, :), ...
        '-',...
        'Color', [1, 0.4, 0.3], ...
        'LineWidth', 2, ...
        'Marker', markers{i}, ...
        'MarkerSize', 10);
end
for i = 1:n2
      % Plot
      h(i,2) = plot(event_time_mat2(i, :),res_mat_2(i, :), ...
        '--',...
        'Color', [0,.6,.6], ...
        'LineWidth', 2, ...
        'Marker', markers{i}, ...
        'MarkerSize', 10);
end


vline(0,'--k');
x_min_plot = -5; % left limit of plot 
xlim([x_min_plot-0.5, x_max+0.5])
ylim([y_min-0.5,y_max+1])
xlabel('time relative to first treated','FontSize',15)
ylabel('treatment effects','FontSize',15)
title('Event Time','FontSize',15)

h1 = ...
legendflex(h(1:n1,1), units_p1, ...
'ref', gca, ...
'anchor', {'nw','nw'}, ...
'buffer', [5 -5], ...
'fontsize',12, ...
'xscale',1.2, ...
'title', 'Quota');
sz = get(h1, 'position');
wid_box = 1.2*sz(3);
sz(3) = wid_box;     % Make the width larger
set(h1,'position', sz);

h2 = ...
legendflex(h(1:n2,2), units_p2, ...
'ref', h1, ...
'anchor', {'sw','nw'}, ...
'buffer', [0 -8], ...
'fontsize', 12, ...
'xscale',1.2, ...
'box','on',...
'title', 'Disclosure');
sz = get(h2, 'position');
sz(3) = wid_box;     % Make the heigth larger
set(h2,'position', sz);
box on
hold off

%% save graphs

saveas(fig,'graph1.png')










