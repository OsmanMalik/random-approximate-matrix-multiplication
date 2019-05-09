function [A, B] = generate_matrices(mat_size, mat_type)
% GENERATE_MATRICES Generate matrix pair (A,B).
%
%   [A, B] = GENERATE_MATRICES(mat_size, mat_type) returns a pair of
%   matrices (A,B), each of size mat_size by mat_size. The argument
%   mat_type is a string that specifies what kind of matrices A and B
%   should be; see code below for a list of all available options.

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

end