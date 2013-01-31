function [obj xstar R M] = dan(X, Y, lambda)
[n p] = size(X);
T = X'*X;
A = [T -T; -T T];
b1 = X'*Y;
b = [b1; -b1];
c = -ones(2*p,1);
[obj sol R M] = parametricsimplex(A,b,c,0,lambda);
R = R(:,1:size(R,2)/2) - R(:,size(R,2)/2+1:end);
xstar = sol(1:p) - sol(p+1:end);
end