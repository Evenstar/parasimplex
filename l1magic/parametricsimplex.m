function [obj xstar R M] = parametricsimplex(A,b,cn,flag, lambda)
R = [];
M = [];
if nargin == 4
    lambda = 0;
end
[m n] = size(A);
A = [A eye(m)];
BB = n+1 : m+n;
NN = 1 : n;
B = A(:,BB);
N = A(:,NN);

xstar = b;
zstar = -cn;
c=[cn;zeros(m,1)];
mu = 1e10;
if flag == 0
    zbar = zeros(size(zstar));
else
    zbar = ones(size(zstar));
end
xbar = ones(size(xstar));
[mu id C] = Max(zstar, zbar, xstar, xbar);
soln = findsoln(n,m,xstar,BB);
R = [R;soln'];
M = [M;mu];
while mu > lambda
    iBN = B\N;
    if C == 1
        j = id;
        ej = circshift(eye(1,size(N,2)),[1,id-1])';
        dx = iBN*ej;
        Xt = dx ./(xstar + mu*xbar);
        Xt(dx == 0) = 0;
        [ux, i] = max(Xt);
        if ux < 1e-9
            disp('unbounded');
            break;
        end
        ei = circshift(eye(1,size(B,2)),[1,i-1])';
        dz = -(iBN)'*ei;
    else
        i = id;
        ei = circshift(eye(1,size(B,2)),[1,id-1])';
        dz = -iBN'*ei;
        Zt = dz./(zstar + mu*zbar);
        Zt(dz == 0) = 0;
        [uz, j] = max(Zt);
        if uz < 1e-9
            disp('unbounded');
            break;
        end
        ej = circshift(eye(1,size(N,2)),[1,j-1])';
        dx = iBN*ej;
    end
    
    t = xstar(i)/dx(i);
    tbar = xbar(i)/dx(i);
    s = zstar(j)/dz(j);
    sbar = zbar(j)/dz(j);
    
    
    
    xstar = xstar - t*dx;
    xbar  = xbar - tbar*dx;
    zstar = zstar - s*dz;
    zbar  = zbar - sbar*dz;
    
    xstar(i) = t;
    xbar(i) = tbar;
    zstar(j) = s;
    zbar(j) = sbar;
    tmp = BB(i);
    BB(i) = NN(j);
    NN(j) = tmp;
    B = A(:,BB);
    N = A(:,NN);
    [mu id C] = Max(zstar, zbar, xstar, xbar);
    modified_xstar = xstar + xbar*mu;
     soln = findsoln(n,m,modified_xstar,BB);
    R = [R;soln'];
    M = [M;mu];

end
obj = c(BB)'*(B\b);

%Only used to obtain solution of parametric simplex
soln = zeros(n,1);
for i = 1 : m
    if BB(i) <= n
        soln(BB(i)) = xstar(i);
    end
end
xstar = soln;
%--Only used to obtain solution of parametric simplex
end

