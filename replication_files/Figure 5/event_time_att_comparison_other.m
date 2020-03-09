
clear
tic 
rng(7)
addpath('functions');
alpha_sig = .05;
% warning('off','all')


%% DATA CLEARNING

data = readtable('data_hr.csv');
% head(data,5)
N = length(unique(table2array(unique(data(:,1))))); % number of units
T1 = length(unique(table2array(unique(data(:,3))))); % number of units
Y = reshape(table2array(data(:,10)),T1,N)'; 
D = reshape(table2array(data(:,4)),T1,N)'; 
D1 = reshape(table2array(data(:,5)),T1,N)'; 
D2 = reshape(table2array(data(:,6)),T1,N)'; 
T = find(sum(D),1)-1; % number of pre-treatment periods
S_max = T1-T; % maximum number of post-treatment periods

S = S_max-23; % truncate data to avoid extrapolating too far
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
id_unit = 7; % unit of interest
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

%% plot1 - absolute level treatment effects for two policies
figure(1)

hold on 
plot([.5,S+.5],[0,0],'--k');
p1 = errorbar(1:S1,gamma_hat_1,ub_1,-lb_1,'Color',...
    [1,0.4,0.3],'LineStyle','none','CapSize',10);
p2 = errorbar(1:S2,gamma_hat_2,ub_2,-lb_2,'Color',...
    [0,.6,.6],'LineStyle','none','CapSize',10);
p3 = plot(1:S1,gamma_hat_1,'Color',[1,0.4,0.3],'LineWidth',2);
p4 = plot(1:S2,gamma_hat_2,'--','Color',[0,.6,.6],'LineWidth',2);

hold off

p1.Marker = 'o';
p1.MarkerSize = 6;
p2.Marker = 'o';
p2.MarkerSize = 6;
p1.LineWidth = 1.5;
p2.LineWidth = 1.5;

xlabel('event time','FontSize',15)
ylabel('ATT estimates','FontSize',15)
xlim([.5 S+.5])
ylim([-10,10])
legend([p3 p4],{'Quota','Disclosure'},'Location',...
    'northwest','FontSize',15);


% %% plot2 - differece 
% figure(2)
% 
% p = errorbar(1:S_min,gamma_hat_diff,ub_vec,-lb_vec);
% 
% hold on
% plot([.5,S+.5],[0,0],'--k')
% hold off
% 
% xlim([.5 12.5])
% ylim([-10,10])
% xlabel('event time')
% ylabel('ATT estimate difference')
% p.Marker = 'o';
% p.LineWidth = 1.5;
% %xlim([.5 S_min+.5])
% lgd = legend('Quota - Disclosure','Location','northwest');
% lgd.FontSize = 12;


%% plot 3 - residual plot
figure(3)

output = att_event_ci(Y,D,S,alpha_sig);
te_mat_hat = output.te_mat_hat;
att_hat = output.att_hat;
ub = output.ub;
lb = output.lb;
res_mat = output.res_mat;
event_time_mat = output.event_time_mat;
B_hat = output.B_hat;

D1_temp = D1;
D2_temp = D2;
ind_temp = sum(D,2)==0;
D1_temp(ind_temp,:) = [];
D2_temp(ind_temp,:) = [];
res_mat_1 = res_mat(~(sum(D1_temp,2)==0),:);
res_mat_2 = res_mat(~(sum(D2_temp,2)==0),:);

% residuals/treatment effects of ever-treated units plot
y_min = min(min(res_mat));
y_max = max(max(res_mat));

subplot(2,1,1)

h = zeros(1,3);
h(1) = plot([0,T+S+1],[0,0],'--k','HandleVisibility','off');


hold on 
p1 = plot(res_mat_1','Color',[1,0.4,0.3],'LineWidth',2);
p2 = plot(res_mat_2','--','Color',[0,.6,.6],'LineWidth',2);
hold off

xlim([1,T+S])
ylim([-20,20])
vline(T,'--k');
xlabel('time','FontSize',15)
ylabel('treatment effects','FontSize',15)
title('Calendar Time','FontSize',15)
lgd = legend([p1(1) p2(1)],'Quota','Disclosure',...
    'Location','northwest','FontSize',12);


% residuals/treatment effects of ever-treated units plot (recentered at
% first treatment time)
subplot(2,1,2)
x_min = min(min(event_time_mat));
x_max = max(max(event_time_mat));
plot([x_min-1,x_max+1],[0,0],'--k')
hold on
if_D1 = ~(sum(D1_temp,2)==0);
for i = 1 : size(res_mat,1)
    if if_D1(i) == 1
        p1 = plot(event_time_mat(i,:),res_mat(i,:),'Color',...
            [1,0.4,0.3],'LineWidth',2);
    else
        p2 = plot(event_time_mat(i,:),res_mat(i,:),'--','Color',...
            [0,.6,.6],'LineWidth',2);
    end
end
hold off
x_min_plot = -5; % left limit of plot 
xlim([x_min_plot, S])
ylim([-20,20])
vline(0,'--k');
xlabel('time relative to first treated','FontSize',15)
ylabel('treatment effects','FontSize',15)
title('Event Time','FontSize',15)
lgd = legend([p1(1) p2(1)],'Quota','Disclosure',...
    'Location','northwest','FontSize',12);


%% save graphs

saveas(figure(1),'graph5a.png')
saveas(figure(3),'graph5b.png')








