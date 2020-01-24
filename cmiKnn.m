%% Conditional Mutual Information
% _cmiKnn_
% Finds the conditional mutual information $I(x;y|z)$ (reads, mutual
% information between $x$ and $y$ given $z$).

%% Citation
%  To cite this code:
%  [ 1 ] Abd AlRahman R. AlMomani, Jie Sun, and Erik Bollt. How Entropic 
%        Regression Beats the Outliers Problem in Nonlinear System 
%        Identification. Chaos 30, 013107 (2020).
%
%
%  [ 2 ] Alexander Kraskov, Harald St¨ogbauer, and Peter Grassberger. 
%        Estimating mutual information. Physical Review E, 69:066-138,2004
%        

%% Inputs
%
% * x : Mxd column wise matrix (source).
% * y : Mxr column wise matrix (distination).
% * z : Mxk column wise matrix of condition set.
% * options: Structure that has two fields: 
%            - options.K: the integer value K for the Kth 
%              nearest neighbor to be used in estimation.
%            - options.distFun: cell array that have name of the 
%              distance function to be used, which should match the
%              built-in distance function defined by matlab (see
%              pdist2)
%

%% Outputs
% _I_ : scalar of the result conditioned mutual information


% The default values for:
%%
%
% * K : 2
% * Distance function: {'minkowski',Inf}


%% Function Body
%
%%
function [I] = cmiKnn(x,y,z)
% Set defaults
distInfo = {'minkowski',Inf}; K = 2;


% If the condition set $z$ is empty, then use the Mutual inforation
% estimator.
if isempty(z)
    [I] = miKSG(x,y);
    return
end
 
% To construct the Joint Space between all variables we have:
JS = cat(2,x,y,z);
% Find the K^th smallest distance in the joint space JS = (x,y,z)
D = pdist2(JS,JS,distInfo{:},'Smallest',K+1)';
epsilon = D(:,end);
% Instead of the above two lines, the one may use the knnsearch function,
% but we found the above implementation is faster.

% Find number of points from $(x,z), (y,z)$, and $(z,z)$ that lies withing the
% K^{th} nearest neighbor distance from each point of themself.
Dxz = pdist2([x,z],[x,z],distInfo{:});
nxz = sum(bsxfun(@lt,Dxz,epsilon),2) - 1;

Dyz = pdist2([y,z],[y,z],distInfo{:});
nyz = sum(bsxfun(@lt,Dyz,epsilon),2) - 1;

Dz = pdist2(z,z,distInfo{:});
nz = sum(bsxfun(@lt,Dz,epsilon),2) - 1;

% VP Estimation formula.
I = psi(K) - mean(psi(nxz+1)+psi(nyz+1)-psi(nz+1)); 