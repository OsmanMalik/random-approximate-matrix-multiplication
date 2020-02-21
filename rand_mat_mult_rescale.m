function C = rand_mat_mult_rescale(A_single, B_single, Y, S, P)
%RAND_MAT_MULT_RESCALE Compute EBC via rescaling Alg. 3 in [Ba16].
%
%C = RAND_MAT_MULT_RESCALE(A_single, B_single, Y, S, P) aims to compute the
%product A_single*B_single by using the EBC defined by the factor matrices
%in Y. Here, S and P should be the sign and permutation functions
%corresponding to a deterministic algorithm. This function uses the 
%"2x O-I" rescaling procedure described in Section 6, Alg. 3, of [Ba16].
%
%REFERENCES:
%
%[Ba16] G. Ballard, A. R. Benson, A. Druinsky, B. Lipschitz, O. Schwartz.
%       Improving the Numerical Stability of Fast Matrix Multiplication.
%       SIAM J. Matrix Anal. Appl. 37(4), pp. 1382-1418, 2016.
        
no_alt  = 2; % Number of outside/inside rescales
N       = size(A_single,1);
Ap      = A_single;
Bp      = B_single;
DA      = eye(N);
DB      = eye(N);

for k = 1:no_alt
    % Outside rescaling step
    DAp     = diag(max(abs(Ap),[],2));
    invDAp  = diag(1./max(abs(Ap),[],2));
    DA      = DA*DAp;
    Ap      = invDAp*Ap;
    
    DBp     = diag(max(abs(Bp),[],1));
    invDBp  = diag(1./max(abs(Bp),[],1));
    DB      = DBp*DB;
    Bp      = Bp*invDBp;
    
    % Inside rescaling step
    D       = diag(sqrt(  max(abs(Bp),[],2).' ./ max(abs(Ap),[],1)  ));
    invD    = diag(sqrt(  max(abs(Ap),[],1) ./ max(abs(Bp),[],2).'  ));
    Ap      = Ap*D;
    Bp      = invD*Bp;
end

% Compute C and scale back
C = rand_mat_mult_C_wrapper(Ap, Bp, Y, 0, S, P);
C = DA*C*DB;

end