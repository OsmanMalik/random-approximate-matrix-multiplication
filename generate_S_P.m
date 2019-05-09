function [S, P] = generate_S_P(max_no_recursions, n, is_random)

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