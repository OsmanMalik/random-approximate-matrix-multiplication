% EXPERIMENT6 Used for computations in Figure 14.

%% Setting
% no_trials: This is the number of products that the average is taken over.
% max_no_recursions: The maximum number of times the algorithm is recursed.
% mat_type: Control what type of matrices are used in the multiplication.
% plot_flag: Control if results are plotted or not.

no_trials = 1e+2;
max_no_recursions = 1; % Not used---leave at 1
mat_type = 'hilbert';
plot_flag = 1;
mat_base_size = 1;

%% Create/load exact and approximate algorithm

Y_EBC = exact_BCRL_decomp((1:7)*1e-1);
[Y_ABC, kappa] = BCRL_decomp(1.5e-2);
n = sqrt(size(Y_EBC{1},1));
mat_size = n^max_no_recursions*mat_base_size;

%% Generate the matrices, compute true C, and do precision conversions

[A, B] = generate_matrices(mat_size, mat_type);
C = A*B;
normC = norm(C, 'fro');
A_single = single(A);
B_single = single(B);

%% Compute C_EBC
[S_det, P_det] = generate_S_P(max_no_recursions, n, false);
C_EBC  = rand_mat_mult_C_wrapper(A_single, B_single, Y_EBC, 0, S_det, P_det);

%% Run computation

% Create matrices that will hold the running average
C_error_fully_random = zeros(max_no_recursions, no_trials);

% Create cell arrays for storing running sum of computed products
C_sum_fully_random = zeros(size(C));

% Create nonrandom S and P
[S_det, P_det] = generate_S_P(max_no_recursions, n, false);

% Main loop
for t = 1:no_trials
    % Create random S and P
    [S_random, P_random] = generate_S_P(max_no_recursions, n, true);
    
    for k = 1:max_no_recursions 
        C_sum_fully_random = C_sum_fully_random + double(rand_mat_mult_C_wrapper(A_single, B_single, Y_ABC, kappa, S_random(1:k), P_random(1:k))); 
        C_error_fully_random(k, t) = norm(C - C_sum_fully_random/t, 'fro')/normC;
    end
end

%% Plot result

if plot_flag == 1
    colors_matlab = get(gca,'colororder');
    plot_colors = colors_matlab(1:3, :);

    % Create plot
    %figure
    p1 = plot([1 no_trials], ones(1,2)*norm(C - double(C_EBC), 'fro')/normC, '--', 'color', 'black', 'linewidth', 2);
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    hold on

    % This is some old code that I currently don't use
    %{
    loglog(C_error_random_S(1,:), 'linewidth', 1, 'linestyle', ':', 'color', plot_colors(1,:))
    loglog(C_error_random_S(2,:), 'linewidth', 1, 'linestyle', ':', 'color', plot_colors(2,:))
    loglog(C_error_random_S(3,:), 'linewidth', 1, 'linestyle', ':', 'color', plot_colors(3,:))

    loglog(C_error_random_P(1,:), 'linewidth', 1, 'linestyle', '-.', 'color', plot_colors(1,:))
    loglog(C_error_random_P(2,:), 'linewidth', 1, 'linestyle', '-.', 'color', plot_colors(2,:))
    loglog(C_error_random_P(3,:), 'linewidth', 1, 'linestyle', '-.', 'color', plot_colors(3,:))
    %}

    p2 = loglog(C_error_fully_random(1,:), 'linewidth', 2, 'linestyle', '-', 'color', plot_colors(1,:));

    ax_min = min(C_error_fully_random(:));
    ax_max = max(C_error_fully_random(:));
    %axis([1 no_trials ax_min*.95 ax_max*1.05])

    legend([p1 p2], {'EBC', 'RandABC'})
    xlabel('Number of trials')
    ylabel('Error of average')

    % Set size of plot
    %x0 = 500;
    %y0 = 500;
    %width = 430;
    %height = 130;
    %set(gcf,'units','points','position',[x0,y0,width,height])
end