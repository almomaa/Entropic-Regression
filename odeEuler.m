function [t,y] = odeEuler(sys, t, y0, h,varargin)

t = (t(1):h:t(2))'; 
n = length(t);
y = repmat(y0,1,n);

for i = 1:n-1
    y(:,i+1) = y(:,i) + h.*sys(t(i),y(:,i),varargin{:});
end
y = y';
end