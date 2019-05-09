function Y_approx = perturb_strassen(Y, noise_level)

Y_approx = Y;
for k = 1:3
    idx = find(Y{k} == 0);
    idx = idx(randsample(length(idx), round(length(idx)/3)));
    idx = [idx; find(Y{k} == 1)];
    Y_approx{k}(idx) = Y_approx{k}(idx) + noise_level*randn(size(idx));
end

end