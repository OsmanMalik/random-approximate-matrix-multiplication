function C = small_mat_mult(A, B, Y)
%SMALL_MAT_MULT Computes A*B using factor matrices in Y
%   
%   C = SMALL_MAT_MULT(A, B, Y) returns an approximation of the product A*B
%   computed using via the bilinear computation defined by the factor
%   matrices in Y, which is a cell containing three factor matrices. Note
%   that A and B must be square n x n matrices, where n is the square root
%   of the number of rows in the factor matrices in Y.
%
%   This code is slow, and not recursive. We do not use it in any of our
%   experiments.

R = size(Y{1},2);
n = sqrt(size(Y{1},1));
C = zeros(n);
for r = 1:R
    q = (Y{1}(:,r).' * A(:)) * (Y{2}(:,r).' * B(:));
    Gamma = reshape(Y{3}(:,r), n, n).';
    C = C + Gamma.*q;
end

end