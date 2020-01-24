function A = polyspace(X, order)

N = size(X,2); A = [];

if (N > 35) || (order <= 4)
    A = highDimension(X,order);
    return;
end

I = getIndexTable(order,N);

for i=1:size(I,1)
    ix = I(i,:); ix(~ix) = [];
    A = cat(2,A,prod(X(:,ix),2));
end

end

function P = getIndexTable( Order,dim )

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


function A = highDimension(X,order)

[m,~] = size(X);
A = ones(m,1);
A = cat(2,A,X); %linear terms

if (order > 1), A = add2ndOrderTerms(A,X); end
if (order > 2), A = add3rdOrderTerms(A,X); end
if (order > 3), A = add4thOrderTerms(A,X); end
if (order > 4), disp('Expensive Request'); end

end

function A = add2ndOrderTerms(A,X)
N = size(X,2);
for i=1:N
    for j=i:N
        A = cat(2,A,prod([X(:,i) X(:,j)],2));
    end
end
end

function A = add3rdOrderTerms(A,X)
N = size(X,2);
for i=1:N
    for j=i:N
        for k=j:N
            A = cat(2,A,prod([X(:,i) X(:,j) X(:,k)],2));
        end
    end
end
end

function A = add4thOrderTerms(A,X)
N = size(X,2);
for i=1:N
    for j=i:N
        for k=j:N
            for u=k:N
                A = cat(2,A,prod([X(:,i) X(:,j) X(:,k) X(:,u)],2));
            end
        end
    end
end
end




