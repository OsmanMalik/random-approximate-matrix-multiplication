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

noise_level = 1e-3;
no_trials = 1e+2;
max_no_recursions = 5;
mat_size = n^max_no_recursions*10;
mat_type = 'normal';

%% Create/load approximate algorithm

strassen_decomp
%laderman_decomp
n = sqrt(size(Y{1},1));

% This is some old code I used when I wanted to try this experiment with an
% actual low-rank 3 by 3 matrix multiplication algorithm we had found.
% I'm not doing that now though since it's extremely slow to use.
%{
load('A_r20_random_search_k_10382')
Y_approx = A;
clear A;
temp1 = 100*max(abs(Y_approx{1}(:)));
temp2 = 100*max(abs(Y_approx{2}(:)));
Y_approx{1} = Y_approx{1}/temp1;
Y_approx{2} = Y_approx{2}/temp2;
Y_approx{3} = temp1*temp2*Y_approx{3};
%}

Y_approx = Y;
for k = 1:3
    idx = find(Y{k} == 0);
    idx = idx(randsample(length(idx), round(length(idx)/3)));
    idx = [idx; find(Y{k} == 1)];
    Y_approx{k}(idx) = Y_approx{k}(idx) + noise_level*randn(size(idx));
end
X = tensor(ktensor(Y)); X = X.data;
X_approx = tensor(ktensor(Y_approx)); X_approx = X_approx.data;
epsilon = sum(X_approx(X==1)-1)/(n^3);

%% Generate the matrices and compute true C

switch mat_type
    case 'hilbert' % Hilbert matrices
        A = hilb(mat_size);
        B = hilb(mat_size);
    case 'binomial' % Binomial matrices
        A = gallery('binomial', mat_size);
        B = A;
    case 'normal' % Standard normal matrices
        A = randn(mat_size);
        B = randn(mat_size);
    case 'normal-normalized' % Normalized standard normal matrices
        A = randn(mat_size);
        A = A/norm(A, 'fro');
        B = randn(mat_size);
        B = B/norm(B, 'fro');
    case 'uniform' % Standard uniform matrices
        A = rand(mat_size);
        B = rand(mat_size);
    case 'adversarial_1' % Type 1 uniform adversarial matrices
        A = rand(mat_size);
        A(:, ceil(mat_size/2):end) = A(:, ceil(mat_size/2):end)/(mat_size^2);
        B = rand(mat_size);
        B(1:floor(mat_size/2), :) = B(1:floor(mat_size/2), :)/(mat_size^2);
    case 'adversarial_2' % Type 2 uniform adversarial matrices
        A = rand(mat_size);
        A(1:floor(mat_size/2), ceil(mat_size/2):end) = A(1:floor(mat_size/2), ceil(mat_size/2):end)*(mat_size^2);
        B = rand(mat_size);
        B(:, 1:floor(mat_size/2)) = B(:, 1:floor(mat_size/2))/(mat_size^2);
    otherwise
        error('Invalid choice for matrix type');
end
C = A*B;
normC = norm(C, 'fro');

%% Run computation

% Create matrix that will store the errors
C_error = zeros(4, max_no_recursions);

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
    C_error(1, k) = norm(C - C_approx_deterministic, 'fro')/normC;
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
        %C_approx_deterministic = rand_mat_mult(A, B, Y_approx, epsilon, S_det(1:k), P_det(1:k), 1);
        %C_approx_random = rand_mat_mult(A, B, Y_approx, epsilon, S_random(1:k), P_random(1:k), 1);
        
        C_approx_fully_random = rand_mat_mult_C_wrapper(A, B, Y_approx, epsilon, S_random(1:k), P_random(1:k));
        %C_approx_random_S = rand_mat_mult_C_wrapper(A, B, Y_approx, epsilon, S_random(1:k), P_det(1:k));
        %C_approx_random_P = rand_mat_mult_C_wrapper(A, B, Y_approx, epsilon, S_det(1:k), P_random(1:k));
        
        C_error(2, k) = C_error(2, k) + norm(C - C_approx_fully_random, 'fro')/(normC*no_trials);
        %C_error(3, k) = C_error(3, k) + norm(C - C_approx_random_S, 'fro')/(normC*no_trials);
        %C_error(4, k) = C_error(4, k) + norm(C - C_approx_random_P, 'fro')/(normC*no_trials);
    end
end

%% Plot results

% Create plot
figure
bar(C_error([1 2], :)')
legend('Deterministic', 'Randomized', 'location', 'northwest')
xlabel('Number of recursions')
ylabel('Error')

% Set size of plot
x0 = 10;
y0 = 10;
width = 430;
height = 130;
set(gcf,'units','points','position',[x0,y0,width,height])