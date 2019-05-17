function [Y_approx, epsilon] = perturb_strassen(Y, noise_level, random_seed)
% PERTURB_STRASSEN Perturb the standard Strassen algorithm.
%
%   [Y_approx, epsilon] = PERTURB_STRASSEN(Y, noise_level, random_seed)
%   takes an exact version of Strassen's algorithm as input Y, with Y being
%   a cell containing the CP factor matrices corresponding to Strassen's
%   algorithm. It then adds mean zero Gaussian noise with standard
%   deviation noise_level to all nonzero elements of the factor matrices.
%   For each factor matrix, it also adds noise to five randomly selected
%   entries that are zero. random_seed is used to set the random seed
%   before adding noise; the random number generator is reshuffled prior to
%   exiting this function. Y_approx is a cell containing
%   the three perturbed factor matrices. epsilon contains the quantity
%   defined in Equation (4) of the paper.

rng(random_seed);

Y_approx = Y;
for k = 1:3
    idx = find(Y{k} == 0);
    idx = idx(randsample(length(idx), round(length(idx)/3)));
    idx = [idx; find(Y{k} == 1)];
    Y_approx{k}(idx) = Y_approx{k}(idx) + noise_level*randn(size(idx));
end

n = sqrt(size(Y{1},1));
X = tensor(ktensor(Y)); X = X.data;
X_approx = tensor(ktensor(Y_approx)); X_approx = X_approx.data;
epsilon = sum(X_approx(X==1)-1)/(n^3);

rng('shuffle')

end