function C = mat_mult_laderman(A, B, nmin)
%MAT_MULT_LADERMAN Compute A*B using Laderman's algorithm [La76]
%
%   C = MAT_MULT_LADERMAN(A, B, nmin) returns the product A*B computed
%   using Laderman's algorithm [La76]. The algorithm uses recursion and can
%   therefore multiply matrices A and B that are of size 3^m. When A and B
%   are of size nmin x nmin or smaller, the matrix multiplication is just
%   done using standard matrix multiply.
%
% REFERENCES:
%   [La76]  J. D. Laderman. A noncommutative algorithm for multiplying
%           3 x 3 matrices using 23 multiplications. Bulletin of the
%           American Mathematical Society, 82(1):126-128, 1976.

N = size(A, 1);

if N <= nmin
    C = A*B;
else
    u = 1 : N/3;
    v = N/3+1 : 2*N/3;
    w = 2*N/3+1 : N;
    
    M1  = mat_mult_laderman(A(u,u) + A(u,v) + A(u,w) - A(v,u) - A(v,v) - A(w,v) - A(w,w), B(v,v), nmin);
    M2  = mat_mult_laderman(A(u,u) - A(v,u), -B(u,v) + B(v,v), nmin);
    M3  = mat_mult_laderman(A(v,v), -B(u,u) + B(u,v) + B(v,u) - B(v,v) - B(v,w) - B(w,u) + B(w,w), nmin);
    M4  = mat_mult_laderman(-A(u,u) + A(v,u) + A(v,v), B(u,u) - B(u,v) + B(v,v), nmin);
    M5  = mat_mult_laderman(A(v,u) + A(v,v), -B(u,u) + B(u,v), nmin);
    M6  = mat_mult_laderman(A(u,u), B(u,u), nmin);
    M7  = mat_mult_laderman(-A(u,u) + A(w,u) + A(w,v), B(u,u) - B(u,w) + B(v,w), nmin);
    M8  = mat_mult_laderman(-A(u,u) + A(w,u), B(u,w) - B(v,w), nmin);
    M9  = mat_mult_laderman(A(w,u) + A(w,v), -B(u,u) + B(u,w), nmin);
    M10 = mat_mult_laderman(A(u,u) + A(u,v) + A(u,w) - A(v,v) - A(v,w) - A(w,u) - A(w,v), B(v,w), nmin);
    M11 = mat_mult_laderman(A(w,v), -B(u,u) + B(u,w) + B(v,u) - B(v,v) - B(v,w) - B(w,u) + B(w,v), nmin);
    M12 = mat_mult_laderman(-A(u,w) + A(w,v) + A(w,w), B(v,v) + B(w,u) - B(w,v), nmin);
    M13 = mat_mult_laderman(A(u,w) - A(w,w), B(v,v) - B(w,v), nmin);
    M14 = mat_mult_laderman(A(u,w), B(w,u), nmin);
    M15 = mat_mult_laderman(A(w,v) + A(w,w), -B(w,u) + B(w,v), nmin);
    M16 = mat_mult_laderman(-A(u,w) + A(v,v) + A(v,w), B(v,w) + B(w,u) - B(w,w), nmin);
    M17 = mat_mult_laderman(A(u,w) - A(v,w), B(v,w) - B(w,w), nmin);
    M18 = mat_mult_laderman(A(v,v) + A(v,w), -B(w,u) + B(w,w), nmin);
    M19 = mat_mult_laderman(A(u,v), B(v,u), nmin);
    M20 = mat_mult_laderman(A(v,w), B(w,v), nmin);
    M21 = mat_mult_laderman(A(v,u), B(u,w), nmin);
    M22 = mat_mult_laderman(A(w,u), B(u,v), nmin);
    M23 = mat_mult_laderman(A(w,w), B(w,w), nmin);
    
    C = [   M6 + M14 + M19 ...
            M1 + M4 + M5 + M6 + M12 + M14 + M15 ...
            M6 + M7 + M9 + M10 + M14 + M16 + M18; ...
            M2 + M3 + M4 + M6 + M14 + M16 + M17 ...
            M2 + M4 + M5 + M6 + M20 ...
            M14 + M16 + M17 + M18 + M21; ...
            M6 + M7 + M8 + M11 + M12 + M13 + M14 ...
            M12 + M13 + M14 + M15 + M22 ...
            M6 + M7 + M8 + M9 + M23                     ];
end

end