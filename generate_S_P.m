function [S, P] = generate_S_P(max_no_recursions, n, is_random)
% GENERATE_S_P Generate random or deterministic sign and permutation
% functions.
%
%   [S, P] = GENERATE_S_P(max_no_recursions, n, is_random) returns sign
%   functions and permutation functions stored in cells. max_no_recursions
%   control how many recursions are considered, i.e., what Q is in our
%   paper. n controls the size of the matrices being multiplied, e.g., 2 in
%   the case of Strassen's algorithm. is_random controls whether or not the
%   returned cells S and P contain random quantities or not: If is_random
%   is false, then the nonrandom sign and permutation functions are
%   returned (i.e., corresponding to deterministic algorithm), and if
%   is_random us true, then random sign and permutation functions are
%   returned.

if is_random
    % Create random S and P
    S = cell(max_no_recursions,1);
    P = cell(max_no_recursions,1);
    for k = 1:max_no_recursions
        S{k} = 2*binornd(1, 0.5, 3, n) - 1;
        P{k} = zeros(3,n);
        for kk = 1:3
            P{k}(kk, :) = randperm(n);
        end
    end
else
    % Create nonrandom S and P
    S = cell(max_no_recursions,1);
    P = cell(max_no_recursions,1);
    for k = 1:max_no_recursions
        S{k} = ones(3, n);
        P{k} = repmat(1:n, 3, 1);
    end
end

end