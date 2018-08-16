%% Strong Basis Selection
% Entropic regression strong basis selection stage (forward step).
% 

%% Version
% This function is a part of Entropic Regression Software Package (erfit),
% version 1.0. To report bugs, comments and suggestions, we appreciate your
% feedback:  Abd AlRahman R. AlMomani, almomaa@clarkson.edu.

%% Function Body
%
%%
function strongIX = sbs( X,y,fwstat, options)

if options.coupFlag
    %If the problem is coupled network, set the node to be 
    %the default influnce source for itself
    strongIX = index;
else
    %Consider the Linear terms 
    strongIX = 2:options.pDim+1;
end %Set default influnce nodes to estimate the tolerence.

Xt = X(:,strongIX);


switch options.sbsMethod
    case 'static'
        I = sort(pcmi(X(:,setdiff(1:size(X,2),strongIX)),y,Xt,options));

        if ~sum(I>0)
            return;
        else
            fwstat.tol = (1-options.alpha)*max(I);
        end

        span = setdiff(1:size(X,2),strongIX);
        for i=span
            if pcmi(X(:,i),y,Xt,options) > fwstat.tol
                Xt = cat(2,Xt,X(:,i));
                strongIX = cat(2,strongIX,i);
            end
        end
        
    case 'dynamic'
        I = pcmi(X(:,setdiff(1:size(X,2),strongIX)),y,Xt,options);

        if ~sum(I>0)
            return;
        else
            fwstat.tol = sCEnt(X,y,Xt,options);%Shuffle test
        end
        pos = strongIX; val = inf; strongIX = [];
        while val > fwstat.tol && ~sum(ismember(pos,strongIX))
            strongIX = cat(2,strongIX,pos);
            I = pcmi(X,y,X(:,strongIX),options);
            [val,pos] = max(I);
            
            fwstat.tol = sCEnt(X(:,pos),y,X(:,strongIX),options);
        end
end
strongIX = sort(strongIX);
assignin('caller','fwstat',fwstat);

%% See Also
% <../html/dataGen.html dataGen>  | <../html/eroptset.html eroptset> |
% <../html/erfit.html erfit> | <../html/cmiVP.html cmiVP>     |
% <../html/miKSG.html miKSG> |  <../html/pcmi.html pcmi> |
% <../html/sbs.html sbs> | <../html/wbr.html wbr> 
% 