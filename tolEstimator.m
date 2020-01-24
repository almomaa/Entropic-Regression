function [tol, I] = tolEstimator( y, options)

N = length(y); I = zeros(1,options.numPerm);

for i=1:options.numPerm
    I(i) = miKnn(y,y(randperm(N)));
end
I = sort(I);
tol = I(ceil(options.alpha*options.numPerm));

