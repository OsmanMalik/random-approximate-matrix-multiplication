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
% mat_size: The size of the matrices multiplied.
% mat_type: Control what type of matrices are used in the multiplication.

no_trials = 1e+2;
max_no_recursions = 5;
mat_size = n^max_no_recursions*10;
mat_type = 'adversarial_1';

%% Create/load exact algorithm

strassen_decomp
%laderman_decomp
n = sqrt(size(Y{1},1));

%% Generate the matrices, compute true C, and do precision conversions

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
A_single = single(A);
B_single = single(B);
C_single = A_single*B_single;

%% Run computation

C_error = zeros(5, max_no_recursions);

% Create nonrandom S and P
S_det = cell(max_no_recursions,1);
P_det = cell(max_no_recursions,1);
for k = 1:max_no_recursions
    S_det{k} = ones(3, n);
    P_det{k} = repmat(1:n, 3, 1);
end

% Compute product using deterministic exact algorithm in single precision
% and compute error
for k = 1:max_no_recursions
    C_approx_deterministic = rand_mat_mult_C_wrapper(A_single, B_single, Y, 0, S_det(1:k), P_det(1:k));
    C_error(1, k) = norm(C - double(C_single), 'fro')/normC;
    C_error(2, k) = norm(C - double(C_approx_deterministic), 'fro')/normC;
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
        
    for k = 1:max_no_recursions  
        C_approx_fully_random = rand_mat_mult_C_wrapper(A_single, B_single, Y, 0, S_random(1:k), P_random(1:k));
        C_approx_random_S = rand_mat_mult_C_wrapper(A_single, B_single, Y, 0, S_random(1:k), P_det(1:k));
        C_approx_random_P = rand_mat_mult_C_wrapper(A_single, B_single, Y, 0, S_det(1:k), P_random(1:k));

        C_error(3, k) = C_error(3, k) + norm(C - double(C_approx_fully_random), 'fro')/(normC*no_trials);
        C_error(4, k) = C_error(4, k) + norm(C - double(C_approx_random_S), 'fro')/(normC*no_trials);
        C_error(5, k) = C_error(5, k) + norm(C - double(C_approx_random_P), 'fro')/(normC*no_trials);
    end
end

%% Plot results

include_refline = false;

figure

bar(C_error(2:end,:)')
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

% Set size
x0 = 10;
y0 = 10;
width = 430;
height = 130;
set(gcf,'units','points','position',[x0,y0,width,height])