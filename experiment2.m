% EXPERIMENT2 Compute and plot error of average for randomized approximate
% matrix product. Used for computations in Figure 1.
%
%   EXPERIMENT2 is a script that computes and plots the error of the
%   average output by the randomized matrix multiplication algorithm.
%   The purpose of this experiment is to verify the claim that our method
%   is correct in expectation for the case of an approximate algorithm.
%   This is done by keeping track of how the error for the average of many
%   independent approximate matrix product changes as the average is taken
%   over more and more approximate matrix products. This is done for
%   different number of recursions.

%% Settings
% noise_level: Control the size of the perturbations added to the exact
%   matrix multiplication algorithm.
% no_trials: This is the number of products that the average is taken over.
% max_no_recursions: The maximum number of times the algorithm is recursed.
% mat_size: The size of the matrices multiplied.
% random_seed: Set the random seed used to perturb Strassen.

noise_level = 1e-3;
no_trials = 1e+2;
max_no_recursions = 3;
mat_type = 'normal';
random_seed = 1;

%% Create/load approximate algorithm

Y = strassen_decomp();
[Y_approx, epsilon] = perturb_strassen(Y, noise_level, random_seed);
n = sqrt(size(Y{1},1));
mat_size = n^max_no_recursions*10;

%% Run computation

% Generate the matrices and compute true C
[A, B] = generate_matrices(mat_size, mat_type);
C = A*B;
normC = norm(C, 'fro');

% Create matrices that will hold the running average
C_error_fully_random = zeros(max_no_recursions, no_trials);
C_error_random_S = zeros(max_no_recursions, no_trials);
C_error_random_P = zeros(max_no_recursions, no_trials);

% Create cell arrays for storing running sum of computed products
C_sum_fully_random = cell(max_no_recursions,1);
C_sum_random_S = cell(max_no_recursions,1);
C_sum_random_P = cell(max_no_recursions,1);
for k = 1:max_no_recursions
    C_sum_fully_random{k} = zeros(size(C));
    C_sum_random_S{k} = zeros(size(C));
    C_sum_random_P{k} = zeros(size(C));
end

% Create nonrandom S and P
S_det = cell(max_no_recursions,1);
P_det = cell(max_no_recursions,1);
for k = 1:max_no_recursions
    S_det{k} = ones(3, n);
    P_det{k} = repmat(1:n, 3, 1);
end

% Main loop
for t = 1:no_trials
    % Create random S and P
    S_random = cell(max_no_recursions,1);
    P_random = cell(max_no_recursions,1);
    for k = 1:max_no_recursions
        S_random{k} = 2*binornd(1, 0.5, 3, n) - 1;
        P_random{k} = zeros(3,n);
        for kk = 1:3
            P_random{k}(kk, :) = randperm(n);
        end
    end
    
    % Compute matrix product and running average
    for k = 1:max_no_recursions
        C_sum_fully_random{k} = C_sum_fully_random{k} + rand_mat_mult_C_wrapper(A, B, Y_approx, epsilon, S_random(1:k), P_random(1:k));
        C_sum_random_S{k} = C_sum_random_S{k} + rand_mat_mult_C_wrapper(A, B, Y_approx, epsilon, S_random(1:k), P_det(1:k));
        C_sum_random_P{k} = C_sum_random_P{k} + rand_mat_mult_C_wrapper(A, B, Y_approx, epsilon, S_det(1:k), P_random(1:k));
        C_error_fully_random(k, t) = norm(C - C_sum_fully_random{k}/t, 'fro')/normC;
        C_error_random_S(k, t) = norm(C - C_sum_random_S{k}/t, 'fro')/normC;
        C_error_random_P(k, t) = norm(C - C_sum_random_P{k}/t, 'fro')/normC;
    end
end

%% Plot results

plot_colors = [62 150 81; 107 76 154; 204 37 41]/255;

% Create plot
figure
p2 = loglog(C_error_fully_random(1,:), 'linewidth', 2, 'linestyle', '-', 'color', plot_colors(1,:));
hold on
p3 = loglog(C_error_fully_random(2,:), 'linewidth', 2, 'linestyle', '-.', 'color', plot_colors(2,:));
p4 = loglog(C_error_fully_random(3,:), 'linewidth', 2, 'linestyle', ':', 'color', plot_colors(3,:));
legend([p2 p3 p4], {'1 recursion', '2 recursions', '3 recursions'})
xlabel('Number of trials')
ylabel('Error of average')

ax_min = min(C_error_fully_random(:));
ax_max = max(C_error_fully_random(:));
axis([1 no_trials ax_min*.95 ax_max*1.05])

% Set size of plot
x0 = 500;
y0 = 500;
width = 430;
height = 130;
set(gcf,'units','points','position',[x0,y0,width,height])