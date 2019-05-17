% EXPERIMENT5 Experiment on impact on error of randomization of exact
% algorithm in single precision. Used for computations in Figures 4 and 5.
%
%   EXPERIMENT5 is a script for running an experiment to verify numerically
%   how the numerical error for an exact algorithm in single precision is
%   impacted by randomization for different number of recursions.

%% Setting
% no_trials: This is the number of products that the average is taken over
%   for the randomized methods.
% max_no_recursions: The maximum number of times the algorithm is recursed.
% mat_type: Control what type of matrices are used in the multiplication.
% plot_flag: Used to control which type of plotting (if any) is used.
% include_refline: If plot_flag is 1, then this controls if a reference
%   line is added or not.

no_trials = 1e+2;
max_no_recursions = 5;
mat_type = 'normal';
plot_flag = 2;
include_refline = true;

%% Create/load exact algorithm

Y = strassen_decomp();
n = sqrt(size(Y{1},1));
mat_size = n^max_no_recursions*10;

%% Generate the matrices, compute true C, and do precision conversions

[A, B] = generate_matrices(mat_size, mat_type);
C = A*B;
normC = norm(C, 'fro');
A_single = single(A);
B_single = single(B);
C_single = A_single*B_single;

%% Run computation

% Create matrix that will store the errors
C_error = nan(5, max_no_recursions, no_trials);

% Create nonrandom S and P
[S_det, P_det] = generate_S_P(max_no_recursions, n, false);

% Compute product using deterministic exact algorithm in single precision
% and compute error
for k = 1:max_no_recursions
    C_approx_deterministic = rand_mat_mult_C_wrapper(A_single, B_single, Y, 0, S_det(1:k), P_det(1:k));
    C_error(1, k, 1) = norm(C - double(C_single), 'fro')/normC;
    C_error(2, k, 1) = norm(C - double(C_approx_deterministic), 'fro')/normC;
end

% Main loop
for t = 1:no_trials    
    % Create random S and P
    [S_random, P_random] = generate_S_P(max_no_recursions, n, true);
        
    for k = 1:max_no_recursions  
        C_approx_fully_random = rand_mat_mult_C_wrapper(A_single, B_single, Y, 0, S_random(1:k), P_random(1:k));
        C_approx_random_S = rand_mat_mult_C_wrapper(A_single, B_single, Y, 0, S_random(1:k), P_det(1:k));
        C_approx_random_P = rand_mat_mult_C_wrapper(A_single, B_single, Y, 0, S_det(1:k), P_random(1:k));

        C_error(3, k, t) = norm(C - double(C_approx_fully_random), 'fro')/normC;
        C_error(4, k, t) = norm(C - double(C_approx_random_S), 'fro')/normC;
        C_error(5, k, t) = norm(C - double(C_approx_random_P), 'fro')/normC;
    end
end

%% Plot results

if plot_flag == 1
    % Compute errors for bar plot
    C_error_plot = zeros(size(C_error, 1), size(C_error, 2));
    C_error_plot(1:2, :) = C_error(1:2, :, 1);
    C_error_plot(3:5, :) = median(C_error(3:5, :, :), 3); %sum(C_error(3:5, :, :), 3)/no_trials;
    C_error_plot = C_error_plot(2:end, :);
    
     % Create plot
    figure
    bar(C_error_plot')
    if include_refline
        hline = refline([0 C_error(1,1)]);
        hline.Color = 'black';
        hline.LineStyle = '--';
        
        legend('Deterministic', 'Fully randomized', 'Random sign', 'Random permutation', 'Standard', 'location', 'northwest')
    else
        legend('Deterministic', 'Fully randomized', 'Random sign', 'Random permutation', 'location', 'northwest')
    end
    xlabel('Number of recursions')
    ylabel('Error')

    % Set size of plot
    x0 = 500;
    y0 = 500;
    width = 430;
    height = 130;
    set(gcf,'units','points','position',[x0,y0,width,height])
elseif plot_flag == 2
    colors_matlab = get(gca,'colororder');
    x_pos = [0 .2 .4 .6];
    bar_width = .15;
    make_boxplots(C_error(2:end, :, :), colors_matlab(1:4, :), {'Deterministic', 'Fully randomized', 'Random sign', 'Random permutation'}, x_pos, bar_width, 'reference_line', C_error(1, 1, 1));
    current_y_lim = get(gca, 'ylim');
    set(gca, 'ylim', [current_y_lim(1)*.1, current_y_lim(2)]);
    
    % Set size of plot
    x0 = 500;
    y0 = 500;
    width = 430;
    height = 130;
    set(gcf,'units','points','position',[x0,y0,width,height])
end