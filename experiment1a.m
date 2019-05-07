% EXPERIMENT1A Test accuracy of Strassen's algorithm for specific random
% configuration.
%
%   EXPERIMENT1A is a script that computes the error in Strassen's
%   algorithm in a setting with limited precision. This is done by
%   multiplying two predetermined 2 by 2 matrices, for a specific sign and
%   permutation configuration. The version of Strassen's algorithm used
%   here is the one described in Section 1.3.11 in [Go13].
%
%   This script utilizes roundsd [Fr19].
%   
% REFERENCES:
%
%   [Go13]  G. H. Golub, C. F. Van Loan. Matrix Computations. 4th Edition.
%           Johns Hopkins University Press, Baltimore, 2013.
%
%   [Fr19]  François Beauducel. Round with significant digits. Available at 
%           Mathwork's website: 
%           https://www.mathworks.com/matlabcentral/fileexchange/26212-round-with-significant-digits?s_tid=srchtitle
%           Accessed 6 May, 2019.

% Add path to roundsd function
addpath('roundsd');

% Create test matrices
A_orig = [.99 .001; .001 .99];
B_orig = [.99 .001; .001 .99];

% Specify S matrices
s1 = [1 1];
s2 = [1 1];
s3 = [1 1];

% Specify P matrices
p1 = [1 2];
p2 = [1 2];
p3 = [1 2];

% Compute matrices M1-M3
I = eye(2);
M1 = I(p1(:),:)*diag(s1);
M2 = I(p2(:),:)*diag(s2);
M3 = I(p3(:),:)*diag(s3);

% Randomized input matrices A, B
A = M1*A_orig*M2.';
B = M2*B_orig*M3.';

% Compute Strassen in limited precision
u = 1; v = 2;

p1h = roundsd((A(u,u)+A(v,v)) * (B(u,u)+B(v,v)), 2);
p2h = roundsd((A(v,u)+A(v,v)) * (B(u,u)), 2);
p3h = roundsd((A(u,u)) * (B(u,v)-B(v,v)), 2);
p4h = roundsd((A(v,v)) * (B(v,u)-B(u,u)), 2);
p5h = roundsd((A(u,u)+A(u,v)) * (B(v,v)), 2);
p6h = roundsd((A(v,u)-A(u,u)) * (B(u,u)+B(u,v)), 2);
p7h = roundsd((A(u,v)-A(v,v)) * (B(v,u)+B(v,v)), 2);

C_inter = [roundsd(p1h+p4h-p5h+p7h, 2) roundsd(p3h+p5h, 2); roundsd(p2h+p4h, 2) roundsd(p1h-p2h+p3h+p6h, 2)]; 

% Compute and output results
C = M1.'*C_inter*M3;
C_error = C - A_orig*B_orig;
fprintf('The error is %.4f\n', norm(C_error, 'fro'));