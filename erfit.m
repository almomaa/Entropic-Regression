%% Entropic Regression
% Entropic regression parameters estimator.
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


%% Syntax
%% 
% 
%   x = erfit(Phi, f)
%   x = erfit(Phi, f, options)
%   [x,Info] = erfit(Phi, f, options)
%   [sol,Info, Mask] = erfit(Phi, f, options)

%% Description
%%
%
% * _x = erfit(Phi, f)_ : Given sampled data for basis functions Phi $\in
% R^{l \times K}$ and   sampled observations,  $f \in R^{l \times d}$,
% erfit finds the sparse solution $x \in R{K \times d}$ that containes the
% true governing parametrs of the dynamic $f$.
% * _x = erfit(Phi, f, options)_ : erfit accept options structure _'options'_
% that controls the computations created by _eroptset_ function.
% * _[x,Info] = erfit(Phi, f, options)_ : _erfit_ provides output structure
% _Info_ that have information about regression process.
% * _[x,Info, Mask] = erfit(Phi, f, options)_ : _erfit_ provides a logical
% matrix _Mask_ where $Mask_{i,j} = 1$ if $x_{i,j} \neq 0$, and $Mask_{i,j}
% = 0$ otherwise.
%

%% Examples
%% 
% Please refer to our main text for the construction of the Phi and f
% matrices. The following is a simplified data construction and function
% call for _erfit_.
%
%   [Phi, f] = dataGen('Lorenz');
%   options = eroptset('sbsMethod','dynamic', 'alpha', 0.95);
%   x = erfit(Phi, f, options);
% 

%% Version
% This function is a part of Entropic Regression Software Package (erfit),
% version 1.0. To report bugs, comments and suggestions, we appreciate your
% feedback:  Abd AlRahman R. AlMomani, almomaa@clarkson.edu.

%% Function Body
%
function [sol,Info, Mask] = erfit(Phi, f, options)

if nargin < 3 %If no options provided, load defaults
    options = eroptset('pDim',size(f,2));
end

%extract system matrcies dimensions
[stat.M,stat.D] = size(Phi); stat.dim = size(f,2);

%Solution placeholder
sol = zeros(stat.D,stat.dim);

fwstat = stat; %Initialize forward stat information.

for i=1:stat.dim
    stat.dimIX = i;
    %Strong Basis Selection (Forward Step).
    strongIX = sbs(Phi, f(:,i), stat, options);
    
    % Update information
    Info(i).fwstat = fwstat;
    Info(i).strIX  = strongIX;

    
    %Weak Basis Removal (Backward Step)
    optimalIX = wbr( Phi(:,strongIX),f(:,i),Info(i).fwstat, options );
    optimalIX = strongIX(optimalIX);
    
    % Update information
    Info(i).bwstat = bwstat;
    Info(i).optIX  = optimalIX;


    
    %If not detected through entropic regression,
    %Include the constant term for influence test.
    if ~ismember(1,optimalIX), optimalIX = cat(2,1,optimalIX); end
    
    %Given the optimal basis, find least squares solution
    sol(optimalIX,i) = nls( Phi(:,optimalIX), f(:,i) );
    
    % Update information
    Info(i).mask = logical(sol(:,i));

end
% Find the solution logical mask required
if nargout == 3, Mask = cat(2,Info.mask); end


%% See Also
% <../html/dataGen.html dataGen>  | <../html/eroptset.html eroptset> |
% <../html/erfit.html erfit> | <../html/cmiVP.html cmiVP>     |
% <../html/miKSG.html miKSG> |  <../html/pcmi.html pcmi> |
% <../html/sbs.html sbs> | <../html/wbr.html wbr> 
% 









