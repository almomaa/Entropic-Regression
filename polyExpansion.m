function [ Aex ,Terms] = polyExpansion( A,Order, timeFlag)
%% Polynomial Expansion 

if nargin < 3
    timeFlag = false;
end

[~,dim] = size(A); Terms = []; Aex = [];

if timeFlag
    if dim <=4
        vec = [' t';' X';' Y';' Z'];
    else
        vec = cat(1,'   t',cat(2,repmat(' X',dim,1),num2str((1:dim)','%02u')));
    end
else
    if dim <=3
        vec = [' X';' Y';' Z'];
    else
        vec = cat(2,repmat(' X',dim,1),num2str((1:dim)','%02u'));
    end
end

P = getIndexTable( Order,dim );
for i=1:size(P,1)
    ix = P(i,:); ix(~ix) = []; str = [];
    Aex = cat(2,Aex,prod(A(:,ix),2));
    
    for j=1:length(ix)
        str = cat(2,str,vec(ix(j),:));
    end
    if isempty(ix)
        str = cat(2,str,'1');
    end

    if dim>9
        Terms = cat(1,Terms,[repmat(' ',1,4*Order-length(str)) str]);
    else
        Terms = cat(1,Terms,[repmat(' ',1,4*Order-length(str)) str]);
    end
end

Terms = mat2cell(Terms,ones(size(Terms,1),1),size(Terms,2));
end

%%
function P = getIndexTable( Order,dim )
%%
ub = (dim+1)^Order - 1; T = [];

for i=0:ub
    T = cat(1,T,dec2base(i, dim+1,Order));
end

T = num2cell(T);
P = cellfun(@encodeChar,T,'UniformOutput',false);
P = cell2mat(P);
G = P(:,2:end)-P(:,1:end-1); 
[ix,~] = find(G<0); P(ix,:) = [];
end

function c = encodeChar(s)
T = 0:100;
B = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ';
ix = strfind(B,s); c = T(ix);
end
