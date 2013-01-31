X = [1 0 1; 2 1 1];
y = [4 9]';
x0 = [0 0 0]';

epsilon = 10;
tic
xp = l1dantzig_pd(x0, X, [], y, epsilon, 5e-2);
toc

%%
epsilon = 10;
lambda = 40;
[n p] = size(X);
G = X'*X;
A = [G, -G, eye(size(G)), zeros(size(G));
     -G, G, zeros(size(G)), eye(size(G))];
B = [eye(p), zeros(p); zeros(p), eye(p)];
N = [G, -G; -G G];
xstar = [X'*y;-X'*y];
zstar = [ ones(p,1); ones(p,1)];
xbar = ones(size(xstar));
zbar = zeros(size(zstar));
Ncal = 1 : 2*p;
Bcal = 2*p+1 : 4*p;
finished = false;

[z1 j] = max(-zstar./zbar.*(zbar>0));
[z2 i] = max(-xstar./xbar.*(xbar>0));
if z1 > z2 
    mu = z1;
    index = Bcal(j);
    cal = 1;
else
    mu = z2;
    index = Ncal(i);
    cal = 2;
end
%while mu > epsilon
    if cal == 1
        ej = zeros(2*p,1); 
        ej(index) = 1;
        deltax = inv(B)*N*ej;
        
        [tmp1 i] = max( deltax ./(xstar + mu*xbar));
        
        ei = zeros(2*p,1);
        ei(Bcal(i)) = 1;
        deltaz = -(inv(B)*N)'*ei;
    else
        ei = zeros(2*p,1);
        ei(index) = 1;
        deltaz = -(inv(B)*N)'*ei;
        
        [tmp2 j] = max( deltaz ./(zstar + mu*zbar));
        
        ej = zeros(2*p,1);
        ej(Ncal(j)) = 1;
        deltax = inv(B)*N*ej;       
    end
    
   t = xstar(Bcal(i))./deltax(Bcal(i));
   tbar = xbar(Bcal(i)) ./deltax(Bcal(i));
   s = zstar(Ncal(j)) ./ deltaz(Ncal(j));
   sbar = zbar(Ncal(j)) ./deltaz(Ncal(j));
   xstar(Ncal(j)) = t;
   xbar(Ncal(j)) = tbar;
   zstar(Bcal(i)) = s;
   zstar(Bcal(j)) = sbar;
   
   xstar = xstar - t*deltax;
   xbar = xbar - tbar*deltax;
   zstar = zstar - s*deltaz;
   zbar = zbar - sbar*deltaz;
%end
