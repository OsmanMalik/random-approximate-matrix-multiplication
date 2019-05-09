function [Y_approx, epsilon] = perturb_strassen(Y, noise_level, random_seed)

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