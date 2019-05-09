% EXPERIMENT3 Experiment on impact on error of randomization of approximate
% algorithm. Used for computations in Figure 2.
%
%   EXPERIMENT3 is a script for running an experiment on how the error
%   changes when randomization is introduced in an approximate bilinear
%   computation for matrix multiplication. The experiment has the following
%   steps:
%
%       1. Create an approximate bilinear computation for matrix
%          multiplication by perturbing Strassen's algorithm.
%       2. Create two test matrices A and B.
%       3. Compute the error for the deterministic algorithm.
%       4. Compute the average error for the randomized algorithm by
%          averaging the results over no_trials trials.
%
%   The error from the deterministic computation is then compared to the
%   average error of the randomized computation in a plot. This experiment
%   is done for multiple different number of recursions.

%% Settings
% noise_level: Control the size of the perturbations added to the exact
%   matrix multiplication algorithm.
% no_trials: This is the number of approximate product computations that
%   the average is taken over. 
% max_no_recursions: The maximum number of times the algorithm is recursed.
% mat_size: The size of the matrices multiplied.
% mat_type: Control what type of matrices are used in the multiplication.
% random_seed: Set the random seed used to perturb Strassen.
% plot_flag: Control if results are plotted or not.

noise_level = 1e-3;
no_trials = 1e+2;
max_no_recursions = 3;
mat_type = 'uniform';
random_seed = 1;
plot_flag = true;

%% Create/load approximate algorithm

Y = strassen_decomp();
[Y_approx, epsilon] = perturb_strassen(Y, noise_level, random_seed);
n = sqrt(size(Y{1},1));
mat_size = n^max_no_recursions*10;

%% Generate the matrices and compute true C

[A, B] = generate_matrices(mat_size, mat_type);
C = A*B;
normC = norm(C, 'fro');

%% Run computation

% Create matrix that will store the errors
C_error = nan(2, max_no_recursions, no_trials);

% Create nonrandom S and P
S_det = cell(max_no_recursions,1);
P_det = cell(max_no_recursions,1);
for k = 1:max_no_recursions
    S_det{k} = ones(3, n);
    P_det{k} = repmat(1:n, 3, 1);
end

% Compute product using deterministic approximate algorithm and compute
% error
for k = 1:max_no_recursions
    C_approx_deterministic = rand_mat_mult_C_wrapper(A, B, Y_approx, epsilon, S_det(1:k), P_det(1:k));
    C_error(1, k, 1) = norm(C - C_approx_deterministic, 'fro')/normC;
end

% Compute product using randomized approximate algorithm no_trial times and
% compute average error
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
    
    for k = 1:max_no_recursions   
        C_approx_fully_random = rand_mat_mult_C_wrapper(A, B, Y_approx, epsilon, S_random(1:k), P_random(1:k));
        C_error(2, k, t) = norm(C - C_approx_fully_random, 'fro')/normC;
    end
end

%% Plot results

if plot_flag
    % Compute errors for bar plot
    C_error_plot = zeros(size(C_error, 1), size(C_error, 2)); 
    C_error_plot(1, :) = C_error(1, :, 1);
    C_error_plot(2, :) = sum(C_error(2, :, :), 3)/no_trials;
    
    % Create plot
    figure
    bar(C_error_plot')
    legend('Deterministic', 'Randomized', 'location', 'northwest')
    xlabel('Number of recursions')
    ylabel('Error')

    % Set size of plot
    x0 = 500;
    y0 = 500;
    width = 430;
    height = 130;
    set(gcf,'units','points','position',[x0,y0,width,height])
end