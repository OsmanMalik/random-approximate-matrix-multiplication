function Y = exact_BCRL_decomp(e_vec)
%exact_BCRL_decomp Create exact CP tensor correspoding to schele in [Bi80].
%
%Y = exact_BCRL_decomp(e_vec) comptues the exact scheme for 12 by 12 matrix
%multiplication discussed in Sec 5 in [Bi80]. This scheme uses a linear
%combination of the inexact APA scheme proposed in [Bi79] to achieve an
%exact scheme. The input e_vec is a vector with errors, which must be of
%length 7 and with distinc elements. These errors are used to construct
%corresponding APA algorithms, which are then combined appropriately into
%an exact scheme.
%
%REFERENCES:
%
%[Bi79] D Bini, M Capovani, F Romani and G Lotti. O(n^{2.7799}) complexity
%       for n x n approximate matrix multiplication. Information Processing
%       Letters 8(5), 1979.
%
%[Bi80] D Bini. Relations between exact and approximate bilinear
%       algorithms. Applications. Calcolo, 1980.

d = 6;
assert(length(e_vec)==d+1, 'e_vec must be of length 7')
e_vec = e_vec(:)'; % Make sure it's a row vector
c = 1000;

V = e_vec.^repelem((0:d)',1,d+1);
b = zeros(d+1,1); 
b(1) = 1;

a = V\b;

Y = cell(3,1);
Y{1} = zeros(12^2, (d+1)*c);
Y{2} = zeros(12^2, (d+1)*c);
Y{3} = zeros(12^2, (d+1)*c);

for k = 1:d+1
    Y_temp = BCRL_decomp(e_vec(k));
    Y{1}(:, 1+(k-1)*c : k*c) = a(k)*Y_temp{1};
    Y{2}(:, 1+(k-1)*c : k*c) = Y_temp{2};
    Y{3}(:, 1+(k-1)*c : k*c) = Y_temp{3};
end

end