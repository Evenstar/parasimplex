%% 
X = [1 1 2 1;
     1 0 1 2;
     0 1 1 0];
Y = [9.1 7.2 2.8]';
lambda = 0.00001;
[obj xstar R M] = dan(X,Y,lambda);
for i = 1 : size(R,1)
    xx = R(i,:)';
    E(i) = max(abs(X'*(Y-X*xx)));
end
addpath('./Optimization');
x0= X'*Y;
xp = l1dantzig_pd(x0, X, [], Y, 0.13, 1e-4)
%%
X = [0 1 1; 1 0.5 0];
Y = [5.2; 2.1];
lambda = 0.000001;
[obj xstar R M] = dan(X,Y,lambda);
for i = 1 : size(R,1)
    xx = R(i,:)';
    E(i) = max(abs(X'*(Y-X*xx)));
end
addpath('./Optimization');
x0= X'*Y;
xp = l1dantzig_pd(x0, X, [], Y, 0.3333, 1e-4)
%%
X = rand(6,10);
beta = randi(5,[10,1]);
Y = X*beta + rand(6,1);
lambda = 1e-5;
[obj xstar R M] = dan(X,Y,lambda)
addpath('./Optimization');
x0= X'*Y;
xp = l1dantzig_pd(x0, X, [], Y, 0.3, 1e-4)
%%
tic
lambda = 0.1;
[obj xstar R M] = dan(X,Y,lambda);
toc
%addpath('./Optimization');
%x0= zeros(size(X,2),1);
%xp = l1dantzig_pd(x0, X, [], Y, lambda, 1e-4)
%%
X = rand(100,200)*3;
beta = randi(10,[200,1]);
Y = rand(100,1)*0.1 + X*beta;
lambda = 0.1;
%%
tic
lambda = 0.1;
[obj xstar R M] = dan(X,Y,lambda);
toc
%%
addpath('./Optimization');
Z=[];
tic
for i = 2 : size(M,1)-1
    i
lambda = M(i);
x0= zeros(size(X,2),1);
xp = l1dantzig_pd(x0, X, [], Y, lambda, 1e-4);
Z=[Z;xp'];
end
toc
%%
lambda = 1014.9
x0= zeros(size(X,2),1);
xp = l1dantzig_pd(x0, X, [], Y, lambda, 1e-4)
%%
X=[2.264 3.678 2.560;
3.678 6.979 4.965;
2.560 4.965 3.724];
[~,~,R,M] = clime(X,3,1e-6);
y = interp1(M,R,0.5)

%%
tic
for i = 1 : 10 : 100
    lambda = R(i);
    xp = l1dantzig_pd(x0, X, [], Y, lambda, 1e-4);
end
toc 
%%
X=[2.264 3.678 2.560;
3.678 6.979 4.965;
2.560 4.965 3.724];
[n p] = size(X);
T = X;
A = [T -T; -T T];
b1 = zeros(n,1); b1(1) = 1;
b = [b1; -b1];
c = -ones(2*p,1);




