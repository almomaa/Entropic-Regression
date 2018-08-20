%% Entropic Regression Options Set
% Create Entropic Regression options structure

%% Syntax
%% 
% 
%   options = eroptset
%   options = eroptset('param1',value1,'param2',value2,...)
% 

%% Description
%%
%
% * options = eroptset (with no input arguments) creates a structure called 
% options that contains the options, or parameters, for the Entropic Regression
% algorithm and sets parameters to [], indicating default values will be used.
% * options = eroptset('param1',value1,'param2',value2,...) creates a structure
% called options and sets the value of 'param1' to value1, 'param2' to 
% value2, and so on. Any unspecified parameters are set to their default values.

%% Options
%
% * _'pDim'_ : Problem dimension: Integer value represent the problem
% dimension. Accept positive integer. Default = 1.
%
% * _'sbsMethod'_ : Strong basis selection method. {'dynamic' or 'static'}.
%
% * _'K"_ : Number of Nearest Neighbors to be used in mutual information
% estimation. Accept positive integer. Default = 2.
%
% * _'distFunc'_ : Distance parametrs to be used in mutual information
% estimation. Accept distance functions supported by Matlab. Default =
% {'minkowski',Inf}.
%
% * _'alpha'_ : Confidence parameter for tolerence etimation through
% shuffle test. alpha in [0,1]. Default = 0.9.
%
% * _'numPerm'_ : Number of permutations to be used in shuffle test. Accept
% positive integer. Default = 100.
%
% * _'coupFlag'_ : Boolean flag. If true, then indicate that the dynamic is
% coupled network of oscillators. Boolean value. Default = false. 
%
% * _'minm'_ : Minimum number of measurements to perform backward step in
% ER. Positive integer. Default = 100.
%

   
%% Examples
%% 
% 
%   options = eroptset
%   options = eroptset('sbsMethod','dynamic')
%   options = eroptset('sbsMethod','dynamic', 'alpha', 0.95)
% 

%% Version
% This function is a part of Entropic Regression Software Package (erfit),
% version 1.0. To report any bugs or issues we appreciate your feedback: 
% Abd AlRahman R. AlMomani, almomaa@clarkson.edu  

%% Function Body
%
%%
function options = eroptset( varargin )
p = inputParser;

% Problem dimention is the dimension of the system before any polynomial
% expansion. Meaning that if the dynamic consists of x,y and z, then pDim =
% 3.
p.addParameter('pDim', 1);

% Strong Basis Selection Method. This parameter decide between two
% approches for selection of best candidate function in the forward step of
% entropic regression. The default method is as described in our main text
% (see citation section in this document) is called 'static' which adopt
% the shuffle test to estimate the tolerence and then it use the same value
% of tolerence through selection process. The second method is the
% 'dynamic' selection, which apply the shuffle test after each successful
% selection of candidate function to apdate the tolerence. Both approches
% have high performance, while the 'static' approch can perform faster than
% the adaptive one specially with high dimensional problems.

p.addParameter('sbsMethod',      'static');

%coupling flag
p.addParameter('coupFlag',  false);



% Alpha
p.addParameter('alpha',   0.90, @(s) (0<=s) && (s<=1));



% numPerm
p.addParameter('numPerm',     10, @isnumeric);

% K - for KSG
p.addParameter('K',   2, @(s) (1<=s));
p.addParameter('distFunc',  {'minkowski',Inf});


p.addParameter('regfun',  @(z,y) z*pinv(z)*y);

p.parse(varargin{:});
options = p.Results;

end

%% See Also
% <../html/dataGen.html dataGen>  | <../html/eroptset.html eroptset> |
% <../html/erfit.html erfit> | <../html/cmiVP.html cmiVP>     |
% <../html/miKSG.html miKSG> |  <../html/pcmi.html pcmi> |
% <../html/sbs.html sbs> | <../html/wbr.html wbr> 
% 
