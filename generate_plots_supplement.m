% GENERATE_PLOTS_SUPPLEMENT Used to generate the plots for supplement.
%
%   GENERATE_PLOTS_SUPPLEMENT is a script we use to generate plots for the
%   supplement. It requires that the results are already stored
%   appropriately in .mat files. This script then simply loads those
%   results, creates the appropriate plots, and save them as FIG files.

results_dir = 'supplement_experiment_results/';
fig_save_dir = 'figures/';
colors_matlab = get(gca,'colororder');
plot_colors = colors_matlab(1:3, :);

%% Make plots for experiment 3

experiment3_files = dir([results_dir, 'experiment3*']);
x_pos = [0, .25];
bar_width = .2;   
for f = 1:length(experiment3_files)
    load([results_dir, experiment3_files(f).name]);
    make_boxplots(C_error, colors_matlab(1:2, :), {'Deterministic', 'Randomized'}, x_pos, bar_width)
    
    % Set size of plot
    x0 = 500;
    y0 = 500;
    width = 430;
    height = 130;
    set(gcf,'units','points','position',[x0,y0,width,height])
    
    savefig([fig_save_dir, 'S-', strrep(experiment3_files(f).name, '_', '-'), '.fig']);
end

%% Make plots for experiment 4

experiment4_files = dir([results_dir, 'experiment4*']);
for f = 1:length(experiment4_files)
    load([results_dir, experiment4_files(f).name]);
    
    % Create plot
    figure
    normC = norm(C, 'fro');
    p1 = plot([1 no_trials], ones(1,2)*norm(C - double(C_single), 'fro')/normC, '--', 'color', 'black', 'linewidth', 2);
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold on

    p2 = loglog(C_error_fully_random(1,:), 'linewidth', 2, 'linestyle', '-', 'color', plot_colors(1,:));
    p3 = loglog(C_error_fully_random(2,:), 'linewidth', 2, 'linestyle', '-.', 'color', plot_colors(2,:));
    p4 = loglog(C_error_fully_random(3,:), 'linewidth', 2, 'linestyle', ':', 'color', plot_colors(3,:));

    ax_min = min(C_error_fully_random(:));
    ax_max = max(C_error_fully_random(:));
    %axis([1 no_trials ax_min*.95 ax_max*1.05])

    legend([p2 p3 p4 p1], {'1 recursion', '2 recursions', '3 recursions', 'Standard'}, 'location', 'west')
    xlabel('Number of trials')
    ylabel('Error of average')

    % Set size of plot
    x0 = 500;
    y0 = 500;
    width = 430;
    height = 130;
    set(gcf,'units','points','position',[x0,y0,width,height])
    
    savefig([fig_save_dir, 'S-', strrep(experiment4_files(f).name, '_', '-'), '.fig']);
end


%% Make plots for experiment 5

experiment5_files = dir([results_dir, 'experiment5*']);
for f = 1:length(experiment5_files)
    load([results_dir, experiment5_files(f).name]);
    
    x_pos = [0 .2 .4 .6 .8];
    bar_width = .15;
    make_boxplots(C_error(2:end, :, :), colors_matlab(1:5, :), {'Deterministic', 'Rescaled 2x O-I', 'Fully randomized', 'Random sign', 'Random permutation'}, x_pos, bar_width, 'reference_line', C_error(1, 1, 1));
    current_y_lim = get(gca, 'ylim');
    set(gca, 'ylim', [current_y_lim(1)*.1, current_y_lim(2)]);
    
    % Set size of plot
    x0 = 500;
    y0 = 500;
    width = 430;
    height = 150;
    set(gcf,'units','points','position',[x0,y0,width,height])
    
    savefig([fig_save_dir, 'S-', strrep(experiment5_files(f).name, '_', '-'), '.fig']);
end
