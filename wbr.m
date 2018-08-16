%% Weak Basis Removal
% Entropic regression weak basis removal stage (backward step).
% 

%% Version
% This function is a part of Entropic Regression Software Package (erfit),
% version 1.0. To report bugs, comments and suggestions, we appreciate your
% feedback:  Abd AlRahman R. AlMomani, almomaa@clarkson.edu.

%% Function Body
%
%%
function optimalIX = wbr( X,y , bwstat, options)

% Enforce high confidance in the backward step
options.alpha = 0.99;
mi = pcmi(X,y,X, options);
if ~sum(mi>0)
    optimalIX = 1:size(X,2);
    assignin('caller','bwstat',bwstat);
    return
else
    mi = [];
end

for i=1:size(X,2)
    mi = cat(1,mi,pcmi(X(:,i),y,X(:,setdiff(1:size(X,2),i)),options));
end

bwstat.tol = (1-options.alpha)*max(mi);

optimalIX = 1:size(X,2);


span = optimalIX;
for i=span
    ix = setdiff(1:size(X,2),i);
    I = pcmi(X(:,i),y,X(:,ix),options);
    if I <= bwstat.tol && sum(logical(optimalIX)) >= 1
        optimalIX(i) = 0;
    end
end
optimalIX = optimalIX(logical(optimalIX));
assignin('caller','bwstat',bwstat);
end

%% See Also
% <../html/dataGen.html dataGen>  | <../html/eroptset.html eroptset> |
% <../html/erfit.html erfit> | <../html/cmiVP.html cmiVP>     |
% <../html/miKSG.html miKSG> |  <../html/pcmi.html pcmi> |
% <../html/sbs.html sbs> | <../html/wbr.html wbr> 
% 