%%
%Just median size. To avoid concentration.
flag = 0;
for i = 1 : 30
m = 10;
n = 20;
A = rand(m,n);
b = rand(m,1);
c = rand(n,1);
lb = zeros(n,1); 
options = optimset('Display','none');
[x fval]= linprog(-c,A,b,[],[],lb,[],[],options);
fval = -fval;
[obj xstar] = parametricsimplex(A,b,c,1);
if abs(fval - obj) > 1e-6
    disp('Wrong');
    flag = 1;
end
% if sum(abs(xstar-x)) > 1e-4
%     disp('Wrong');
%     flag = 1;
% end
end
if flag == 0
disp('Tests passed');
else 
    disp('Tests Failed');
end
cvx_begin
   variable x(n)
   dual variables y z
   minimize( -c' * x  )
   subject to
      y : A * x <= b;
      
      z : x >= 0;
cvx_end
-obj
%%
for i = 1 : 30
m = 10;
n = 20;
A = randi(100,[m,n]);
b = randi(100,[m,1]);
c = randi(200,[n,1]);
lb = zeros(n,1);
options = optimset('Display','none');
[x fval]= linprog(-c,A,b,[],[],lb,[],[],options);
fval = -fval;
[obj xstar] = parametricsimplex(A,b,c,1);
if abs(fval - obj) > 1e-6
    disp('Wrong');
end
if sum(abs(xstar-x)) > 1e-4
    disp('Wrong');
end
end
disp('Tests passed');

%%
[obj xstar ] = parametricsimplex(A,b,c,1);
options = optimset('Display','none');
[x fval]= linprog(-c,A,b,[],[],lb,[],[],options);

%%



for i = 1 : 30
m = 10;
n = 20;
A = randi(100,[m,n])-5;
b = randi(100,[m,1])-5;
c = randi(200,[n,1])-10;
lb = zeros(n,1);
options = optimset('Display','none');
[x fval]= linprog(-c,A,b,[],[],lb,[],[],options);
fval = -fval;
[obj xstar] = parametricsimplex(A,b,c,1);
if abs(fval - obj) > 1e-6
    disp('Wrong');
end
if sum(abs(xstar-x)) > 1e-4
    disp('Wrong');
end
end
disp('Done');
%%
%this test shows that specifically set lambda, the 
%pd interior point method produces the same result as
%ps.
X = [1 1 2 1;
     1 0 1 2;
     0 1 1 0];
Y = [9.1 7.2 2.8]';
%X = rand(10,[40,500]);
%Y = rand(10,[40,1]);
%X = rand(4,5);
%Y = rand(4,1);
[n p] = size(X);
T = X'*X;
A = [T -T; -T T];
b1 = X'*Y;
lambda = 0.1587;
b = [lambda+b1; lambda-b1];
c = -ones(2*p,1);
tic
[obj sol R M] = parametricsimplex(A,b,c,0,1e-6);
toc
R = R(:,1:size(R,2)/2) - R(:,size(R,2)/2+1:end);
R(end,:)

addpath('./Optimization');
tic
x0= X'*Y;
xp = l1dantzig_pd(x0, X, [], Y, lambda, 1e-6)
toc



%%
%cvx
m = 100000;
n = 3;
A = rand(m,n);
b = rand(m,1);
c = rand(n,1);
tic
cvx_begin
   variable x(n)
   dual variables y z
   minimize( -c' * x  )
   subject to
      y : A * x <= b;     
      z : x >= 0;
cvx_end
toc

tic
[x fval]= linprog(-c,A,b,[],[],lb,[],[],options);
toc
















