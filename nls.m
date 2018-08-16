function x = nls( A, y )

% This function simplify the solution of the problem y = Ax
% if A and y given to be the problem with true basis stored
% in A, the nls return the Least Squares solution.
% If x is solution optained by a method such as compressed sensing
% where we have the solution have very low parameters but are not equal 
% to zero... then
% x = nls( x_cs )
% will return the solution in sparse form by thresholding the parameters
% based on the energy ... where any parameter has energy less than 0.1% of
% the mean energy of the parameters will be thresholded.
% | x_cs_i |^2  < 0.001 * mean(x_cs.^2)
% This function is useful to compare different methods solution in term of
% sparse structure.
if nargin < 2
    x = A;
else
    x = pinv(A)*y;
end
P = zeros(size(x));

for i=1:length(x)
    P(i) = abs(x(i)).^2/mean(x.^2);
end

ix = find(P<=1e-3); x(ix) = 0;

end
