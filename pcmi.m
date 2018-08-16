%% Pairwise Conditional Mutual Information
% The  _pcmi_ function solves the equation:
%
% $$ I_{i} = C_{X_{i} \rightarrow y | Z} $$
%
% _vec = pcmi(X, y, Z, options)_ 
% Finds the pairwise conditional mutual information between each column in
% $X$ and $y$ conditioned on the columns of $Z$. This function replace the
% traditional _for_ loop through the dynamic reccurenc and returns a row
% vector _vec_ with the results.
% For more information in causation entropy, please refer to:
%  [ 1 ] rio erik 2014 siam


%% Inputs
%%
%
% * X : Mxd column wise matrix (source).
% * y : Mx1 column vector (distination).
% * Z : Mxk column wise matrix of condition set.
%

%% Outputs
% _vec_ : 1xd row vector of the result conditioned mutual information

%% Options
% _options_ is structure contain options to pass to the mutual information
% estimator. The default for _options_ structure is to be created by
% _eroptset_ function. The required field at this version is _options.K_
% which the number of nearest neighbors to be used in the mutual
% information estimator.
%

%% Function Body
%
%
function vec = pcmi(X, y, Z, options)

if size(X,2)==1
    vec = cmiVP(X,y,feval(options.regfun,Z,y), options.K);
else
    vec = cat(2,pcmi(X(:,1:end-1), y, Z, options),...
                cmiVP(X(:,end),y,feval(options.regfun,Z,y), options.K));
                
end

%% See Also
% <../html/dataGen.html dataGen>  | <../html/eroptset.html eroptset> |
% <../html/erfit.html erfit> | <../html/cmiVP.html cmiVP>     |
% <../html/miKSG.html miKSG> |  <../html/pcmi.html pcmi> |
% <../html/sbs.html sbs> | <../html/wbr.html wbr> 
%

