%% KSG Mutual Information Estimator
%

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
%% 
% _miKnn_
% Finds the mutual information $I(x;y)$ 
% (reads, mutual information between $x$ and $y$).


%% Inputs
%%
%
% * x : Mxd column wise matrix (source).
% * y : Mxr column wise matrix (distination).
% * options: Structure that has two fields: 
%            - options.K: the integer value K for the Kth 
%              nearest neighbor to be used in estimation.
%            - options.distFun: cell array that have name of the 
%              distance function to be used, which should match the
%              built-in distance function defined by matlab (see
%              pdist2)
%

%% Outputs
% _I_ : scalar of the result estimated mutual information



%% Examples
%%
%
% * I = miKSG(x,y)               ; finds the mutual information I(x;y)
% using default values.
% * options.K = 3; I = miKSG(x,y,options)  ; finds I(x;y) using the 3rd NN.
% * options.distFun = {'euclidean'}; 
%   I = miKSG(x,y,options); 
%   finds I(x;y) using euclidean distance
%
% The default values for:
%%
%
% * K : 2
% * Distance function: {'minkowski',Inf}


%% Function Body
%
%%
function I = miKnn(x,y)
% Initialize and verify inputs.
distInfo = {'minkowski',Inf}; K = 2;


% To construct the Joint Space between all variables we have:
JS = [x,y];  n = size(JS,1);


% Find the K^th smallest distance in the joint space JS
D = pdist2(JS,JS,distInfo{:},'Smallest',K+1)';
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
I = psi(K) + psi(n) - mean(psi(nx+1)+psi(ny+1));