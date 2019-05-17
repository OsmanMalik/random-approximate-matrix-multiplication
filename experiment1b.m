% EXPERIMENT1B Compute average accuracy of randomized Strassen's algorithm.
% Used for computations in Section 2.2 of paper.
%
%   EXPERIMENT1B is a script that computes the average, or expected value,
%   of a randomized version of Strassen's algorithm in a setting with
%   limited precision. This is done by multiplying two predetermined 2 by 2
%   matrices. We can compute the expected value exactly since there are
%   only 64*8 possible values for the random variables to take, and they're
%   all equiprobable. The version of Strassen's algorithm used here is the
%   one described in Section 1.3.11 in [Go13]. 
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

% Construct test matrices
A_orig = [.99 .001; .001 .99];
B_orig = [.99 .001; .001 .99];

% Create cells to store all possible sign and permutation configurations
S = cell(64,1);
P = cell(8,1);

% Create all 64 sign configurations
for k1 = 1:2
for k2 = 1:2
for k3 = 1:2
for k4 = 1:2
for k5 = 1:2
for k6 = 1:2
    idx = k6 + 2*(k5-1) + 2^2*(k4-1) + 2^3*(k3-1) + 2^4*(k2-1) + 2^5*(k1-1);
    S{idx} = [(-1)^k1 (-1)^k2; (-1)^k3 (-1)^k4; (-1)^k5 (-1)^k6];
end
end
end
end
end
end

% Create all 8 permutation configurations
P{1} = [1 2; 1 2; 1 2];
P{2} = [2 1; 1 2; 1 2];
P{3} = [1 2; 2 1; 1 2];
P{4} = [2 1; 2 1; 1 2];
P{5} = [1 2; 1 2; 2 1];
P{6} = [2 1; 1 2; 2 1];
P{7} = [1 2; 2 1; 2 1];
P{8} = [2 1; 2 1; 2 1];

% Create matrix that will contain the product A*B.
C_sum = zeros(size(A_orig));

% Compute the product
for Sidx = 1:64
    for Pidx = 1:8
        % Construct M matrices
        s1 = S{Sidx}(1,:);
        s2 = S{Sidx}(2,:);
        s3 = S{Sidx}(3,:);
        p1 = P{Pidx}(1,:);
        p2 = P{Pidx}(2,:);
        p3 = P{Pidx}(3,:);
        M1 = I(p1(:),:)*diag(s1);
        M2 = I(p2(:),:)*diag(s2);
        M3 = I(p3(:),:)*diag(s3);
        
        % Compute randomized versions of matrices
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
        C_inter = roundsd([p1h+p4h-p5h+p7h p3h+p5h; p2h+p4h p1h-p2h+p3h+p6h], 2);
        
        % "Derandomize" the computed matrix product
        C_rounded = M1.'*C_inter*M3;
        
        % Add up all the computed products
        C_sum = C_sum + C_rounded;
    end
end

% Compute error of expectation compared to true (double precision) product
C_avg = C_sum/(64*8);
C_error = C_avg - A_orig*B_orig;
fprintf('The error is %.4f\n', norm(C_error, 'fro'));