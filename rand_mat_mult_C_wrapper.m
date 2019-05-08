function C = rand_mat_mult_C_wrapper(A, B, Y, epsilon, S, P)
% RAND_MAT_MULT_C_WRAPPER Computes A*B approximately using either
% rand_mat_mult_c_single or rand_mat_mult_c_double depending on the inputs
% A and B.
%
%   C = RAND_MAT_MULT_C_WRAPPER(A, B, Y, epsilon, S, P) returns an
%   approximation of the product A*B computed via a bilinear computation
%   defined by the factor matrices in Y, the adjustment factor epsilon, the
%   realization of Rademacher random variables in S, and the realization of
%   permutation random variables in P. If either A or B are single, then A,
%   B, Y and epsilon are all converted to single; otherwise, they are all
%   converted to double. The entries of the S and P are always converted to
%   int32.
%
%   Compile the source code in rand_mat_mult_C.c using the script
%   compile_rand_mat_mult_C.m to get the properly named functions that are
%   required for this function.

%% Choose appropriate C function

if isa(A, 'single') || isa(B, 'single')
    A = single(A);
    B = single(B);
    for k = 1:length(Y)
        Y{k} = single(Y{k});
    end
    epsilon = single(epsilon);
    for k = 1:length(S)
        S{k} = int32(S{k});
        P{k} = int32(P{k});
    end
    C = rand_mat_mult_C_single(A, B, Y, epsilon, S, P);
else
    A = double(A);
    B = double(B);
    for k = 1:length(Y)
        Y{k} = double(Y{k});
    end
    epsilon = double(epsilon);
    for k = 1:length(S)
        S{k} = int32(S{k});
        P{k} = int32(P{k});
    end
    C = rand_mat_mult_C_double(A, B, Y, epsilon, S, P);
end

end