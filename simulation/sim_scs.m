
clear
tic 
rng(7)
restoredefaultpath
addpath('functions','data');
alpha_sig = .05;
% warning('off','all')



%% LOAD DATA
data = load('sim_data.mat');

%% PREPARATION
case_names = fieldnames(data);
ATT_struct = struct();
ATT_event_struct = struct();

%% START SIMULATION
for f = 1:numel(case_names)

    case_name = case_names{f};
    data_case = data.(case_name);   % each data_case contains one N-T combination
    num_sims = size(data_case, 3);     % # of simulation
    fprintf('Processing case: %s, nsims = %d\n', case_name, num_sims);

    % Prepare storage for this case
    ATT_event = cell(num_sims,1);
    
    %%PARALLEL LOOP OVER SIMULATIONS
    parfor i = 1:num_sims
        t0 = tic

        % Print progress every 200 simulations
        if mod(i,200) == 0
            disp(['  Simulation ', num2str(num_sims-i), '/', num2str(num_sims), ' completed...']);
        end

        %%Extract simulation i
        data_i = data_case(:,:,i);
        
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
        output_t = att_event_ci(Y,D,S,alpha_sig);
        att_hat_t = output_t.att_hat;
        
        %%CONFIDENCE INTERVALS
        ci_l_t = att_hat_t - output_t.ub;
        ci_u_t = att_hat_t - output_t.lb;
        
        %%COLLECT RESULTS
        ATT_event{i} = [(1:S)', att_hat_t - eff(1:S), att_hat_t, (att_hat_t-eff(1:S)).^2, (ci_l_t<=eff(1:S)) & (ci_u_t>=eff(1:S))];
        
        %%USED TIME PER SIM
        elapsed_time(i,1) = toc(t0);  
        fprintf('Sim %d | elapsed %.2f sec\n', ...
            i,  elapsed_time(i,1));
    end
    
    %%STORE RESULTS IN STRUCT
    ATT_event_struct.(case_name) = ATT_event;
end
disp('Finished')
%%
save('sim_result_SSC.mat', 'ATT_struct', 'ATT_event_struct')

%% SUMMARY
%%extract results
sim_mat = fieldnames(ATT_event_struct);
sim_mat = sim_mat{1};
sim_mat = ATT_event_struct.(sim_mat)(cellfun(@(x) size(x,1) == 10, ATT_event_struct.(sim_mat)));
sim_mat = cat(3, sim_mat{:});

%%summary
ATT = squeeze(sim_mat(:, 3, :));
ATT_Var = var(ATT, 0, 2); 


elapsed_time = repmat(sum(elapsed_time), length(ATT_Var), 1);

summary = mean(sim_mat, 3);
summary(:,3) = ATT_Var;    
summary = [summary, elapsed_time];
summary = num2cell(summary);
summary = cell2table(summary, 'VariableNames', ...
    {'event_time','bias','var', 'mse','coverage_rate', 'time_total(s)'});
%% SAVE
if ~exist('output', 'dir')
    mkdir('output')
end
writetable(summary, 'output/sim_summary_SSC.csv');






