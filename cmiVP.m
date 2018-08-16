%% Conditional Mutual Information
% _cmiVP_
% Finds the conditional mutual information $I(x;y|z)$ (reads, mutual
% information between $x$ and $y$ given $z$).
% The _cmiVP_ estimator adopt the principle of nearest neighbor... more

%% Citation
%
%  To cite this code:
%  [ 1 ] Abd AlRahman R. AlMomani, Jie Sun, and Erik Bollt. How Entropic 
%        Regression Beats the Outliers Problem in Nonlinear System 
%        Identification. ... Preprint.
%

%% Inputs
%
% * x : Mxd column wise matrix (source).
% * y : Mxr column wise matrix (distination).
% * z : Mxk column wise matrix of condition set.
%
% Optional Inputs:
% _cmiVP_ accept any built-in distance function defined by matlab in
% addition to the user defined distance functions. please see
% <matlab:pdist2 pdist2> documentation for more info.
%
% _cmiVP_ accept also the integer value K for the Kth nearest neighbor to
% be used in estimation.

%% Outputs
% _I_est_ : scalar of the result conditioned mutual information



%% Examples
%%
%
% * I = cmiVP(x,y,z)               ; finds the mutual information I(x;y|z)
% using default values.
% * I = cmiVP(x,y,z,3)             ; finds I(x;y|z) using the 3rd NN.
% * I = cmiVP(x,y,z,'cityblock')   ; finds I(x;y|z) using cityblock distance.
% * I = cmiVP(x,y,z,10,'euclidean'); finds I(x;y|z) using euclidean distance
% and the 10th NN.
%
% The default values for:
%%
%
% * K : 2
% * Distance function: {'minkowski',Inf}


%% Function Body
%
%%
function [I_est] = cmiVP(x,y,z,varargin)

distInfo = {'minkowski',Inf}; K = 2;
if nargin > 3
    if isnumeric(varargin{1})
        K = varargin{1};
        
    elseif strcmp(varargin{1},'seuclidean') || ...
            strcmp(varargin{1},'mahalanobis') || ...
            strcmp(varargin{1},'minkowski')
        
        distInfo = {varargin{1},varargin{2}};
        if nargin > 5
            K = varargin{3};
        end
    else
       distInfo = varargin(1);
       if nargin > 4
           K = varargin{2};
       end
    end
end

% If the condition set $z$ is empty, then use the Mutual inforation
% estimator.
if isempty(z)
    [I_est] = miKSG(x,y,options);
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
I_est = psi(K) - mean(psi(nxz+1)+psi(nyz+1)-psi(nz+1)); 



%% See Also
% <../html/dataGen.html dataGen>  | <../html/eroptset.html eroptset> |
% <../html/erfit.html erfit> | <../html/cmiVP.html cmiVP>     |
% <../html/miKSG.html miKSG> |  <../html/pcmi.html pcmi> |
% <../html/sbs.html sbs> | <../html/wbr.html wbr> 
% 