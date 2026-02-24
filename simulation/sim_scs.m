%% SET ENVIRONMENT
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clear
rng(7)
restoredefaultpath
addpath('functions','data');
% warning('off','all')



%% LOAD DATA
data = load('sim_data.mat');


%% START SIMULATION
t0 = tic;

%%preparation
case_names = fieldnames(data);
ATT_event_struct = struct();

%%estimation
for f = 1:numel(case_names)

    case_name = case_names{f};
    data_case = data.(case_name);   % each data_case contains one N-T-r combination
    num_sims = size(data_case, 3);     % # of simulation
    fprintf('Processing case: %s, nsims = %d\n', case_name, num_sims);

    % Prepare storage for this case
    ATT_event = cell(num_sims,1);
    
    %%PARALLEL LOOP OVER SIMULATIONS
    parfor i = 1:num_sims
        t0_i = tic


        %%Extract simulation i
        data_i = data_case(:,:,i);
        
        eff = (1:10)';  % true event-tim effect 
        att = unique(data_i(:,5));  % true ATT
        %%RESHAPE DATA
        N = length(unique(data_i(:,1)));
        T1 = length(unique(data_i(:,2)));
        Y = reshape(data_i(:,3),T1,N)'; 
        D = reshape(data_i(:,4),T1,N)'; 
        
        %%DEFINE TREATMENT MATRIX
        T = find(sum(D,1),1)-1;
        S_max = T1-T;
        S = S_max;
        T1 = T + S;
        Y = Y(:,1:T1);
        D = D(:,1:T1);
        D_S = D(:,T+1:T+S);
        
        %%ESTIMATION         
        output_t = att_event(Y,D,S);
        elapsed_time(i,1) = toc(t0_i);
        
        
        %%COLLECT RESULTS
        att_hat_t = output_t.att_hat;
        ATT_event{i} = att_hat_t;
        
        %%USED TIME PER SIM  
        fprintf('Sim %d | elapsed %.2f sec\n', ...
            i,  elapsed_time(i,1));
    end
    
    %%STORE RESULTS IN STRUCT
    ATT_event_struct.(case_name) = ATT_event;

    %%Intermidiate save after each case for safety
    % save('sim_result_SSC.mat', 'ATT_struct', 'ATT_event_struct')
end
disp('Finished')
toc(t0)


%% EXTRACT RESULT AND SAVE

%%create folder or delete existing results
if ~exist('output', 'dir')
    mkdir('output')
else
    delete(fullfile('output','*ssc*'))
end


%%extract and save
%%one csv file for each case, one simulation per row, event time ATT per column
for f = 1:numel(case_names)
    case_name = case_names{f};

    %%extract results
    sim_result = ATT_event_struct.(case_name); 
    sim_result = cat(2, sim_result{:}).';
    varNames = strcat("event_time_", string(1: width(sim_result)));
    sim_result = [sim_result, elapsed_time];
    varNames = [varNames, "time_sec"];
    sim_result = array2table(sim_result, 'VariableNames', varNames);

    %%save
    file_name = ['output/sim_results_', case_name, '_ssc.csv'];
    writetable(sim_result, file_name);
    fprintf('Case %s | saved in %s\n', case_name, file_name);
end







