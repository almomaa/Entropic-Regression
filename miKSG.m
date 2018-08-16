%% KSG Mutual Information Estimator
%

%% Copyright and Citation
%  Copyright 2018, Entropic Regression Package
%
% <https://webspace.clarkson.edu/~ebollt/Website-C3S2/index.html Clarkson Center for Complex System Science (C3S2)> 
%
%
%  To cite this code:
%  [ 1 ] Abd AlRahman R. AlMomani, Jie Sun, and Erik Bollt. How Entropic 
%        Regression Beats the Outliers Problem in Nonlinear System 
%        Identification. ... Preprint
%
%
%  [ 2 ] Alexander Kraskov, Harald St¨ogbauer, and Peter Grassberger. Estimating 
%        mutual information. Physical Review E - Statistical Physics, Plasmas,
%        Fluids, and Related Interdisciplinary Topics, 2004.


%% 
% _miKSG_
% Finds the conditional mutual information $I(x;y)$ (reads, mutual
% information between $x$ and $y$).
% The _miKSG_ estimator adopt the principle of nearest neighbor... more



%% Inputs
%%
%
% * x : Mxd column wise matrix (source).
% * y : Mxr column wise matrix (distination).
%
%
% Optional Inputs:
% _miKSG_ accept any built-in distance function defined by matlab in
% addition to the user defined distance functions. please see
% <matlab:pdist2 pdist2> documentation for more info.
%
% _miKSG_ accept also the integer value K for the Kth nearest neighbor to
% be used in estimation.

%% Outputs
% _I_ : scalar of the result estimated mutual information



%% Examples
%%
%
% * I = miKSG(x,y)               ; finds the mutual information I(x;y)
% using default values.
% * I = miKSG(x,y,3)             ; finds I(x;y) using the 3rd NN.
% * I = miKSG(x,y,'cityblock')   ; finds I(x;y) using cityblock distance.
% * I = miKSG(x,y,10,'euclidean'); finds I(x;y) using euclidean distance
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
function I = miKSG(x,y,varargin)
% Initialize and verify inputs.
distInfo = {'minkowski',Inf}; k = 2;
if nargin > 2
    if isnumeric(varargin{1})
        k = varargin{1};
        
    elseif strcmp(varargin{1},'seuclidean') || ...
            strcmp(varargin{1},'mahalanobis') || ...
            strcmp(varargin{1},'minkowski')
        
        distInfo = {varargin{1},varargin{2}};
        if nargin > 4
            k = varargin{3};
        end
    else
       distInfo = varargin(1);
       if nargin > 3
           k = varargin{2};
       end
    end
end


% To construct the Joint Space between all variables we have:
JS = [x,y];  n = size(JS,1);


% Find the K^th smallest distance in the joint space JS
D = pdist2(JS,JS,distInfo{:},'Smallest',k+1)';
epsilon = D(:,end); %Set threshold value

% Find points on x with pairwise distance
% less than threshold value
Dx = pdist2(x,x,distInfo{:});
nx = sum(bsxfun(@lt,Dx,epsilon),2) - 1;

% Find points on y with pairwise distance
% less than threshold value
Dy = pdist2(y,y,distInfo{:});
ny = sum(bsxfun(@lt,Dy,epsilon),2) - 1;

% KSG Estimation formula.
I = psi(k) + psi(n) - mean(psi(nx+1)+psi(ny+1));


%% See Also
% <../html/dataGen.html dataGen>  | <../html/eroptset.html eroptset> |
% <../html/erfit.html erfit> | <../html/cmiVP.html cmiVP>     |
% <../html/miKSG.html miKSG> |  <../html/pcmi.html pcmi> |
% <../html/sbs.html sbs> | <../html/wbr.html wbr> 
%