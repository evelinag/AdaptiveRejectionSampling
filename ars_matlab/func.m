function ys = func(xs)
    
compositeWeights = [0.25; 0.25; 0.25; 0.25];
pis = [0.25; 0.25; 0.25; 0.25];
a0 = 2.0;
b0 = 20.0;

ys = zeros(size(xs));
for ix=1:length(xs)
    x = xs(ix);

    ys(ix) = (a0 - 1.0) * log(x) - b0 * x ...
        + gammaln(sum(compositeWeights) * x) ...
        - sum(gammaln(x * compositeWeights)) ...
        + sum((x * compositeWeights - 1) .* (log(pis)));
end
