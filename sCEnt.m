%% Shuffle test with Causation Entropy
% Shuffle test.
% 

%% Version
% This function is a part of Entropic Regression Software Package (erfit),
% version 1.0. To report bugs, comments and suggestions, we appreciate your
% feedback:  Abd AlRahman R. AlMomani, almomaa@clarkson.edu.

%% Function Body
%
%%
function tol = sCEnt( X, y, z, options )

%If no conditioning set provided, perform shuffle test 
%with transfer entropy
if isempty(z), tol = sTEnt(X,y, options); return; end

[T,~] = size(y); tol = [];
for i=1:options.numPerm
    ix = randperm(T);
    tol = cat(1,tol,pcmi(X, y(ix),z,options));
end

tol = sort(tol,1);
tol = sort(tol(ceil(options.alpha*size(tol,1)),:));
tol = tol(ceil(options.alpha*length(tol)));
end