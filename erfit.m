%% Entropic Regression
%
% 

%% Citation
%  To cite this code:
%  [ 1 ] Abd AlRahman R. AlMomani, Jie Sun, and Erik Bollt. How Entropic 
%        Regression Beats the Outliers Problem in Nonlinear System 
%        Identification. Chaos 30, 013107 (2020).
%        
%% Description
% S = erfit(Phi, Xdot): solve the inverse problem:
% given Xdot and Phi, find S such that Xdot = Phi*S.
% erfit returns a matrix S of coefficient for a linear regression of 
% the responses in vector field Xdot as a linear combination of the 
% basis functions in the matrix Phi. 
% 


%% Inputs
%%
%
% * Phi  : Nxd basis matrix (i.e. polynomial expansion of the state variable).
% * Xdot : Nxr vector field (derivative of the state variables).
% * options: Structure that has the erfit options. See eroptset.m
%

%% Outputs
% S : dxr coefficients matrix.

%% 
% Author: Abd AlRahman R. AlMomani
% Clarkson University, 2019
% For any comments and feedback please email the author at:
% aaalmoma@clarkson.edu, or, almomaniar@gmail.com

%% Function Body
%
%%
function S = erfit(Phi, Xdot)

% Tolerence estimation paraters:
% The algorithm is not 'very' sensitive to these parameters
% and they can be neglected in case of using fixed tolerence
% as discussed below.
options.alpha   = 0.99;
options.numPerm = 500;


% Least squares fitting (svd, minimum energy)
lsfit  = @(ix,iy) pinv(Phi(:,ix))*Xdot(:,iy);

dim = size(Xdot,2); [~,N] = size(Phi);

%Initialize the solution
S = zeros(N,dim);

% Start Entropic Regression
for i=1:dim
    % default tolerence initialization
    tol = tolEstimator(Xdot(:,i),options);
    % Note: the above function may be the most expensive step
    % and in case of large number of measurement (>2000),
    % we may replace the above function by fixed tol such as
    % tol = 0.01. 
    
    % For the ith dimension:
    IX = 1:N; %Explore all the candidate functions
    
    % Select the strong candidate functions from IX
    % through the forward ER.
    IX = erForward(Phi,Xdot(:,i),IX,tol);
    
    % Eliminate the weak candidate functions from IX
    % through the backward ER.
    IX = erBackward(Phi, Xdot(:,i), IX, tol);

    % Compute the Least squares solution for the selected functions.
    S(IX,i) = lsfit(IX,i);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Forward Entropic Regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IX = erForward(A,y,IXin, tol)
val = inf; ix = []; IX = [];

N = size(A,2);

while val>tol %Start the forward ER
    IX = cat(2,IX,IXin(ix)); %Add the selected strong candidate. 
    D = -inf(1,N);      %Initialize mutual information vector

    for i=1:N
        if ~ismember(IXin(i),IX) %Don't recheck what already selected
            % Find the extra information added by the ith candidate 
            % given the strong candidates selected so far
            D(i) = cmiKnn(A(:,[IX IXin(i)])*pinv(A(:,[IX IXin(i)]))*y,...
                         y,...
                         A(:,IX)*pinv(A(:,IX))*y);
        end
    end

    % Find the strongest candidate that have the maximum extra information
    [val,ix] = max(D);

    %If the maximum extra information (val) is more than the minimum
    %accepted influnce (tol)... then the function will be added to the
    %strong candidates (see IX = [IX,IXin(ix)]; above ). Otherwise,
    %terminate the search.
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Backward Entropic Regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IX = erBackward(A, y, IX, tol)
val = -inf; ix = [];

while ( val < tol ) %Start the backward ER
    IX(ix) = [];           %Eliminate the selected weak candidate
    D = inf(1,length(IX)); %Initialize mutual information vector
    for i=1:length(D)      %For the ith strong candidate
        rem = setdiff(IX,IX(i)); %find all other candidates except i
        % Find the extra information added by the ith strong candidate 
        % given the all other strong candidates. (causation entropy). 
        D(i) = cmiKnn(A(:,IX)*pinv(A(:,IX))*y,y,...
                     A(:,rem)*pinv(A(:,rem))*y);
    end
    % Find the weakest candidate that have the minimum extra information
    [val,ix] = min(D);
    
    %If the minimum extra information (val) is less than the minimum
    %accepted influnce (tol)... then the function will be removed from the
    %strong candidates (see IX(ix) = []; above ), because that means 
    %the function is too weak and has no significant influnce compared to
    %the other candidates. Otherwise, that if val>tol, that mean the 
    %weakest candidate is strong, and have significant influnce, So,
    %terminate the backward step.
end
end