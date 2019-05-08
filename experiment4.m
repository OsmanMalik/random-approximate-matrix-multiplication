% EXPERIMENT4 Compute and plot error of average for randomized exact single
% precision matrix product.
%
%   EXPERIMENT4 is a script to run an experiment with the purpose of is to
%   exploring what the expectation of our randomized algorithms is in the
%   setting of an exact fast algorithm but with numerical rounding. This is
%   done by keeping track of how the error for the average of many
%   independent matrix products changes as the average is taken over more
%   and more approximate matrix products. This is done for different number
%   of recursions.

%% Setting
% no_trials: This is the number of products that the average is taken over.
% max_no_recursions: The maximum number of times the algorithm is recursed.
% mat_size: The size of the matrices multiplied.
% mat_type: Control what type of matrices are used in the multiplication.

no_trials = 1e+4;
max_no_recursions = 5;
mat_size = n^max_no_recursions*10;
mat_type = 'normal';

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
    
    for k = 1:max_no_recursions
        %C_sum{k} = C_sum{k} + double(rand_mat_mult(A_single, B_single, Y, 0, S_random(1:k), P_random(1:k), 1));
        
        C_sum_fully_random{k} = C_sum_fully_random{k} + double(rand_mat_mult_C_wrapper(A_single, B_single, Y, 0, S_random(1:k), P_random(1:k))); 
        %C_sum_random_S{k} = C_sum_random_S{k} + double(rand_mat_mult_C_wrapper(A_single, B_single, Y, 0, S_random(1:k), P_det(1:k)));
        %C_sum_random_P{k} = C_sum_random_P{k} + double(rand_mat_mult_C_wrapper(A_single, B_single, Y, 0, S_det(1:k), P_random(1:k)));
        
        C_error_fully_random(k, t) = norm(C - C_sum_fully_random{k}/t, 'fro');
        %C_error_random_S(k, t) = norm(C - C_sum_random_S{k}/t, 'fro');
        %C_error_random_P(k, t) = norm(C - C_sum_random_P{k}/t, 'fro');
    end
end

%% Plot result

plot_colors = [62 150 81; 107 76 154; 204 37 41]/255;

% Create plot
figure
p1 = plot([1 no_trials], ones(1,2)*norm(C - C_single, 'fro'), '--', 'color', 'black');
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
p3 = loglog(C_error_fully_random(2,:), 'linewidth', 2, 'linestyle', '-.', 'color', plot_colors(2,:));
p4 = loglog(C_error_fully_random(3,:), 'linewidth', 2, 'linestyle', ':', 'color', plot_colors(3,:));

ax_min = min(C_error_fully_random(:));
ax_max = max(C_error_fully_random(:));
%axis([1 no_trials ax_min*.95 ax_max*1.05])

legend([p1 p2 p3 p4], {'Standard', '1 recursion', '2 recursions', '3 recursions'})
xlabel('Number of trials')
ylabel('Error of average')

% Set size of plot
x0 = 10;
y0 = 10;
width = 430;
height = 130;
set(gcf,'units','points','position',[x0,y0,width,height])