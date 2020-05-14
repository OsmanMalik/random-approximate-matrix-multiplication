function [Y, varargout] = BCRL_decomp(e)
%BCRL_decomp Create the CP tensor corresponding to scheme in [Bi79].
%
%Y = BCRL_decomp(e) computes the scheme for 12 by 12 matrix multiplication
%discussed in [Bi79] and Sec 5 of [Bi80]. The output Y is on the same
%format as strassen_decomp and laderman_decomp, so it can be used with
%rand_mat_mult_C_wrapper.
%
%As far as we can see, there are some typos in Eq (5.2) in [Bi80]; we have
%fixed these in the code below. Specifically, there are minus signs missing
%in three places in the definintion of W322 in that Eq in [Bi80].
%
%REFERENCES:
%
%[Bi79] D Bini, M Capovani, F Romani and G Lotti. O(n^{2.7799}) complexity
%       for n x n approximate matrix multiplication. Information Processing
%       Letters 8(5), 1979.
%
%[Bi80] D Bini. Relations between exact and approximate bilinear
%       algorithms. Applications. Calcolo, 1980.

f = 1/e;
tkr3 = @(a,b,c) tkron(tkron(a,b),c);
vec = @(x) x(:);

% [U,V,W] matrices for (3,2,2) scheme

U322 = ...
    [1 0 1 0 1 0 0 0 0 0
     0 0 0 0 0 1 1 0 1 0 
     0 0 0 0 0 0 0 0 e e
     0 0 0 e e 0 0 0 0 0
     1 1 0 1 0 0 0 0 0 0 
     0 0 0 0 0 1 0 1 0 1];

V322 = ...
    [e  0  0 -e  0  1 -1  1  0  1
     0  0  0  0  e  0 -1  0  1  0
     0 -1  0  1  0  0  0  0  0  e
     1 -1  1  0  1  e  0  0 -e  0];
V322 = V322([1 3 2 4], :); % To avoid having to transpose B
 
W322 = ... % Note: W322 in [BCRL] is incorrect.
    [f  f -f  f  0  0  0  0  0  0
     0  0 -f  0  f  0  0  0  0  0
     0  0  0  1  0  1  0  0  0 -1
     1  0  0  0 -1  0  0  0  1  0
     0  0  0  0  0  0  0 -f  0  f
     0  0  0  0  0  f  f -f  f  0];
W322 = W322([1 3 5 2 4 6],:); % To avoid having to transpose C 

% [U,V,W] matrices for (2,2,3) scheme

U223 = V322([1 3 2 4], :);
V223 = U322([1 4 2 5 3 6], :);
W223 = W322([1 4 2 5 3 6], :);

% [U,V,W] matrices for (2,3,2) scheme

U232 = U322([1 4 2 5 3 6], :);
V232 = W322;
W232 = V322;

% Combine into [U,V,W] matrices for (12,12,12) scheme
U = matricize(tkr3(tensorize(U322,3,2), tensorize(U232,2,3), tensorize(U223,2,2)));
V = matricize(tkr3(tensorize(V322,2,2), tensorize(V232,3,2), tensorize(V223,2,3)));
W = matricize(tkr3(tensorize(W322,3,2), tensorize(W232,2,2), tensorize(W223,2,3)));

% Correct W so that it is of the same format as output from strassen_decomp
% and laderman_decomp
M = reshape(1:12^2,12,12);
id(vec(M')) = 1:length(vec(M'));
W = W(id,:);

% Put together return cell
Y = cell(3,1);
Y{1} = U;
Y{2} = V;
Y{3} = W;

% Compute epsilon (kappa) when necessary
if nargout == 2
    n = 12;
    X = tensor(ktensor(Y));
    epsilon = 0;
    for i = 1:n
        for j = 1:n
            for l = 1:n
                epsilon = epsilon + X(i+n*(l-1), l+n*(j-1), j+n*(i-1)) - 1;
            end
        end
    end
    epsilon = epsilon/(n^3);
    varargout{1} = epsilon;
end

end

%% Some help functions

function T = tensorize(A, m, n)
T = reshape(A,m,n,size(A,2));
end

function A = matricize(T)
A = reshape(T,size(T,1)*size(T,2),size(T,3));
end

function C = tkron(A,B)

szA = num2cell(size(A));
szB = num2cell(size(B));
C = repelem(A,szB{:}) .* repmat(B,szA{:});

end
