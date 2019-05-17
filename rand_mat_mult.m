function C = rand_mat_mult(A, B, Y, epsilon, S, P, rec_level)
%RAND_MAT_MULT Computes A*B approximately
%
%   C = RAND_MAT_MULT(A, B, Y, epsilon, S, P, rec_level) returns an
%   approximation of the product A*B computed via a bilinear computation
%   defined by the factor matrices in Y, the adjustment factor epsilon, the
%   realization of Rademacher random variables in S, the realization of 
%   permutation random variables in P, and the current recursive level
%   rec_level.
%
%   This code is very slow compared, and we therefore do not use it in any
%   experiments. The C code is fast and is used instead.

%% Initialize C to single or double precision array depending on input

if isa(A, 'single') || isa(B, 'single')
    C = single(zeros(size(A)));
else
    C = zeros(size(A));
end

%% Compute product using recursion

N = size(A,1);
R = size(Y{1},2);
n = sqrt(size(Y{1},1));

if rec_level > length(S) % recursive level exceeds available random variables; use standard matrix multiply
    C = A*B;
else % run randomized scheme on this recursive level
    s = S{rec_level};
    p = P{rec_level};
    for i = 1:n
        for j = 1:n
            C_block = zeros(N/n);
            for r = 1:R
                if Y{3}((p(1,i)-1)*n+p(3,j), r) ~= 0
                    A_block = zeros(N/n);
                    B_block = zeros(N/n);
                    for k = 1:n
                        for L = 1:n
                            if Y{1}(p(1,k)+(p(2,L)-1)*n, r) ~= 0
                                A_block = A_block + Y{1}(p(1,k)+(p(2,L)-1)*n, r) * s(1,k) * s(2,L) * A(1+(k-1)*N/n:k*N/n, 1+(L-1)*N/n:L*N/n);
                            end
                            if Y{2}(p(2,k)+(p(3,L)-1)*n, r) ~= 0
                                B_block = B_block + Y{2}(p(2,k)+(p(3,L)-1)*n, r) * s(2,k) * s(3,L) * B(1+(k-1)*N/n:k*N/n, 1+(L-1)*N/n:L*N/n);
                            end
                        end
                    end
                    C_block = C_block + Y{3}((p(1,i)-1)*n+p(3,j), r) * rand_mat_mult(A_block, B_block, Y, epsilon, S, P, rec_level+1);
                end
            end
            C(1+(i-1)*N/n:i*N/n, 1+(j-1)*N/n:j*N/n) = 1/(1+epsilon) * s(1,i) * s(3,j) * C_block;
        end
    end
end

end