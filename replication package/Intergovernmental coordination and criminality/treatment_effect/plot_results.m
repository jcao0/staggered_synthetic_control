% This Matlab script plot Figure for the empirical replication of:
% "Increasing Intergovernmental Coordination to Fight Crime: Evidence from Mexico"

% The Figure compare results estimated by staggered synthetic control and
% general synthetic control (Xu, 2017)

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clear
t0 = tic;
rng(7)
restoredefaultpath
addpath('output');
alpha_sig = .05;



%% load results
results_ssc = readtable("results_ssc.csv");
results_gsc = readtable("results_gsc.csv");
outcomes = {'hom_all_rate','hom_ym_rate','theft_violent_rate','theft_nonviolent_rate','presence_strength','co_num','war'};

%% plot figure
labels = ['a':'z']; 
titles = {'Homicide rate', 'Cartel-related homicide rate','Violent theft rate','Nonviolent theft rate',  'Cartel strength', 'Number of cartels', 'Cartel war'};
f = figure('Units','pixels','Position',[100,100,800,1000]);
nCols = 2;
nRows = ceil(length(outcomes) / nCols);
for i =1:length(outcomes)
    outcome = string(outcomes{i});
    results_gsc_i = results_gsc(results_gsc.outcome==outcome,:);
    S_need = max(results_gsc_i.event_time);
    results_ssc_i = results_ssc(results_ssc.outcome==outcome & results_ssc.eventTime<=S_need,:);
    
    gamma_hat_ssc = results_ssc_i.attEstimate;
    ub_ssc = results_ssc_i.attEstimate - results_ssc_i.confidenceInterval_l;
    lb_ssc = results_ssc_i.attEstimate - results_ssc_i.confidenceInterval_u;

    gamma_hat_gsc = results_gsc_i.ATT;
    ub_gsc = results_gsc_i.ATT - results_gsc_i.CI_lower;
    lb_gsc = results_gsc_i.ATT - results_gsc_i.CI_upper;

    subplot(nRows, nCols, i)     % 3 rows, 2 columns, k-th plot
    hold on
    plot([-.5,S_need-.5],[0,0],'--k');


    x = 0:S_need-1; % set event-time to start from zero
    fill([x fliplr(x)], ...
        [gamma_hat_ssc-ub_ssc; flipud(gamma_hat_ssc-lb_ssc)], ...
        [1,0.4,0.3], 'FaceAlpha', 0.1, 'EdgeColor','none');

    fill([x fliplr(x)], [gamma_hat_gsc-ub_gsc; flipud(gamma_hat_gsc-lb_gsc)], ...
        [0,0.6,0.6], 'FaceAlpha', 0.1, 'EdgeColor','none');

    p3 = plot(x,gamma_hat_ssc,'Color',[1,0.4,0.3],'LineWidth',2);
    p4 = plot(x,gamma_hat_gsc,'--','Color',[0,.6,.6],'LineWidth',2);

    % title and label
    title = replace(titles{i}, "_", " ");
    label = strcat('(', labels(i), ") ",title);
    
    ax = gca;
    pos = ax.Position;  % [left bottom width height] in normalized figure units
    annotation('textbox', [pos(1), pos(2)-0.08, pos(3), 0.05], ...
           'String', label, ...
           'HorizontalAlignment','center', ...
           'VerticalAlignment','middle', ...
           'LineStyle','none', 'FontSize',14)
    xlabel('event time','FontSize',12)
    ylabel('ATT estimates','FontSize',12)
    xlim([-.5 S_need-.5])
    legend([p3 p4],{'SSC','GSC'},'Location',...
        'northwest','FontSize',10);
    
    % Reduce height slightly to create more vertical space
    pos = get(gca,'Position');   % current axes position, [left bottom width height]
    pos(4) = pos(4)*0.9;  % make axes 90% of gsc height
    set(gca,'Position',pos)
    hold off
end

%% save figure
if ~exist('output', 'dir')
    mkdir('output')
end
print(f, 'output/Figure3_application_results.png', '-dpng', '-r300')  % high-resolution PNG

%%save running time
elapsed_seconds = toc(t0);
output_file = 'output/running_time_plot_results.txt';
fid = fopen(output_file, 'w');
fprintf(fid, '"Running Time (seconds)"\n');
fprintf(fid, '%.6f\n', elapsed_seconds);
fclose(fid);
