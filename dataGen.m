%% Data Generator
% Data Generator. Generates data (dynamic and basis expansion) for chosen
% dynamical system.
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
%        Identification. ... Preprint.
%


%% Syntax
%% 
% 
%   [Phi, f] = dataGen( system )
%   [Phi, f] = dataGen( system, 'param1',value1,'param2',value2,...'paramN',valueN )
%   [Phi, f, Info] = dataGen( system, 'param1',value1,... )
%

%% Description
%%
%
% * _[Phi, f] = dataGen( system )_ : Generates the fifth order polynomial
% expansion (Phi) and the dynamic (f) for the system given by the
% charachter vector (system).
%
% * _[Phi, f] = dataGen( system, 'param1',value1, 'param2',value2, ...,
% 'paramN',valueN )_ : Generates the basis expansion (Phi) and the dynamic
% (f) for the system given by the charachter vector (system), and sets the
% data generation control parameters 'param1' to value1, 'param2' to
% value2, and so on. Any unspecified parameters are set to their default
% values.
%
% * _[Phi, f, Info] = dataGen( system, 'param1',value1,... )_ : Return a
% structure Info that contains detailed information about the data
% generation process.
%

%% Options (input pairs)
%
% * _'system'_ : is the only "required" input which is a string array with
% the name of the system to be generated. system takes one of the values:
% {'Lorenz'}. See the Note is version section.
%
% * _'NoiseLevel'_: The standard deviation of the noise added to the
% system. Refere to our main text for $\epsilon_1$. _'NoiseLevel'_ is a
% pair input should be followed with double value. The default value = 0.
%
% * _'CorruptionStd'_: The standard deviation of the corruption noise added
% to the system. Refere to our main text for $\epsilon_2$.
% _'CorruptionStd'_ is a pair input should be followed with double value.
% The default value = 0.
%
% * _'CorruptedProb'_: The standard deviation of the corruption noise added
% to the system. Refere to our main text for $p$. _'CorruptedProb'_ is a
% pair input should be followed with double value. The default value = 0. 
%
% * _'SampleSize'_: Number of sample points (data) to be generated,
% excluding the first 1000 points which will be considered as transient
% points. _'SampleSize'_ is a pair input should be followed with integer
% value. The default value = 500.  
%
% * _'ExpansionOrder'_: The expansion orde used to construct the basis
% functions matrix Phi. _'ExpansionOrder'_ is a pair input should be
% followed with integer  value. The default value = 5. (Warning: dataGen
% function is prepaired to take expansion order up to 35, but with high
% dimentional systems the construction of such high dimension matrix can be
% slow, or even out of memory size).
%
% * _'initCond'_ : The initial condition for the system to be generated. If
% ignored, the default is to generate the initial condition randomly inside
% the attractor.
%
% * _'stepSize'_ : The step size to be used with odeEuler integrator. 
%

%% Examples
%% 
% Please refer to our main text for the construction of the Phi and f
% matrices. The following is a simplified data construction and function
% call for _dataGen_.
%
%   [Phi, f] = dataGen('Lorenz'); 
%   [Phi, f, Info] = dataGen('Lorenz','SampleSize',1000,'NoiseLevel',0.05);
%   [Phi, f] = dataGen('Lorenz','SampleSize',1000,...
%                               'NoiseLevel',0.05,...
%                               'CorruptedProb',0.1,...
%                               'CorruptedStd',10);
% 

%% Version
% This function is a part of Entropic Regression Software Package (erfit),
% version 1.0. To report bugs, comments and suggestions, we appreciate your
% feedback:  Abd AlRahman R. AlMomani, almomaa@clarkson.edu.
%
% Note: In this primary version of dataGen function, we introduce the
% basics of data generation process through Lorenz system example. We are
% working on preparing large datasets number to be released in future
% versions. Your suggestion is highly appreciated. 
%

%% Function Body
%
function [A,Y, Info] = dataGen( system , varargin)

% Inputs verifications and initialization
p = inputParser;
p.addRequired('system',@isstr);

p.addParameter('NoiseLevel',      0, @isnumeric);
p.addParameter('SampleSize',    500, @isnumeric);
p.addParameter('ExpansionOrder',  5, @isnumeric);
p.addParameter('CorruptedProb',   0, @(s) (0<=s) && (s<=1));
p.addParameter('CorruptionStd',   0, @isnumeric);
p.addParameter('initCond',        0, @isnumeric);
p.addParameter('stepSize',        0, @isnumeric);

p.parse(system,varargin{:}); inputs = p.Results;

system   = inputs.system;          %System name
Noise    = inputs.NoiseLevel;      %Noise STD
T        = inputs.SampleSize;      %Sample Size
exOrder  = inputs.ExpansionOrder;  %Expansion Order
CorRatio = inputs.CorruptedProb;   %Corruption Probability
CorrStd  = inputs.CorruptionStd;   %Corruption STD
corrFlag = (0<CorRatio);           %Flag indicates if data is corrupted
x0       = inputs.initCond;        %initial condition.
h        = inputs.stepSize;        %step size.

switch system
                
    case 'Lorenz'
        if ~x0, x0 = [1 1 1]' + 0.1*randn(3,1); end
        if ~h , h  = 0.015;      end
        
        sig = 10; beta = 8/3; rho = 28;
        
        %can accept N initial conditions
        lorenzODE = @(t,x) [(sig.*(x(2:3:end)-x(1:3:end)));
                            (x(1:3:end).*(rho-x(3:3:end))-x(2:3:end));
                            (x(1:3:end).*x(2:3:end)-beta.*x(3:3:end))];
        %skip 1000 transients                
        tspan = [0 T*h+1000*h];
        [~,Z] = odeEuler(lorenzODE, tspan, x0,h);

        Z = Z(1001:end,:);
        X = Z(1:end-1,:);
        Y = (1/h).*(Z(2:end,:)-Z(1:end-1,:));
        
        P = zeros(7,size(Y,2));
        P(1:7,:) =  [ 0    0     0;
                    -10   28     0;
                     10   -1     0;
                      0    0  -8/3;
                      0    0     0;
                      0    0     1;
                      0   -1     0];
                  
                  
end


[A , Terms] = polyExpansion( X,exOrder);
P = cat(1,P,zeros(size(A,2)-size(P,1),size(Y,2)));

if corrFlag
    corrCount = ceil(CorRatio*T);
    cix       = randperm(T,corrCount); %corrpted index
    rix       = 1:T; %all index
else
    corrCount = 0; cix = []; rix = 1:T;
end

normalized = @(q) q;
Y(cix,:) = Y(cix,:) + CorrStd.*normalized(randn(size(Y(cix,:))));
Y(rix,:) = Y(rix,:) + Noise.*normalized(randn(size(Y(rix,:))));

if nargout > 2
    Info.x0         = x0;
    Info.h          = h;
    Info.P          = P;
    Info.X          = X;
    Info.Y          = Y;
    Info.A          = A;
    Info.Terms      = Terms;
    Info.cix        = cix;
    Info.rix        = rix;
    Info.corrCount  = corrCount;
end

%% See Also
% <../html/dataGen.html dataGen>  | <../html/eroptset.html eroptset> |
% <../html/erfit.html erfit> | <../html/cmiVP.html cmiVP>     |
% <../html/miKSG.html miKSG> |  <../html/pcmi.html pcmi> |
% <../html/sbs.html sbs> | <../html/wbr.html wbr> 
% 
