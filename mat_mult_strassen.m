function C = mat_mult_strassen(A, B, nmin)
%MAT_MULT_STRASSEN  Computes A*B using Strassen's algorithm
%
%   C = MAT_MULT_STRASSEN(A, B, nmin) returns the product A*B computed
%   using Strassen's algorithm, specifically Algorithm 1.3.1 in [1].
%   The algorithm uses recursion and can therefore multiply matrices A and
%   B that are of size 2^m. When A and B are of size nmin x nmin or
%   smaller, the matrix multiplication is just done using standard matrix
%   multiply.
%
% REFERENCES:
%   [1] G. H. Golub, and C. F. Van Loan. Matrix Computations. Fourth 
%       edition, The Johns Hopkins University Press, 2013.

N = size(A, 1);

if N <= nmin
    C = A*B;
else
    u = 1:N/2;
    v = N/2+1:N;

    M1 = mat_mult_strassen(A(u,u)+A(v,v), B(u,u)+B(v,v), nmin);
    M2 = mat_mult_strassen(A(v,u)+A(v,v), B(u,u), nmin);
    M3 = mat_mult_strassen(A(u,u), B(u,v)-B(v,v), nmin);
    M4 = mat_mult_strassen(A(v,v), B(v,u)-B(u,u), nmin);
    M5 = mat_mult_strassen(A(u,u)+A(u,v), B(v,v), nmin);
    M6 = mat_mult_strassen(A(v,u)-A(u,u), B(u,u)+B(u,v), nmin);
    M7 = mat_mult_strassen(A(u,v)-A(v,v), B(v,u)+B(v,v), nmin);

    C = [M1+M4-M5+M7 M3+M5; M2+M4 M1-M2+M3+M6]; 
end

end