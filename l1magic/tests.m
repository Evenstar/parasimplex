%Example of dantzig selector
load('diabetes.mat');
obj = dantzig(x,y);
figure(1);
X = sum(abs(obj.M),2);
for i = 1 : 10
plot(X, obj.M(:,i));hold all;
end
title('dantzig selector for diabetes dataset');
ylabel('Coefficients');
xlabel('|\beta|_1');

%%
%Example of clime estimator, covariance is tridiagonal,with diagonal all
%ones, and sub-diagonal all 0.5s. Trivial case where n > p.
clear all
m = 5;
N = 50;
S=diag(ones(2*m+1,1)) + 0.5*diag(ones(2*m,1),1) + 0.5*diag(ones(2*m,1),-1);
X = mvnrnd(zeros(1,2*m+1), S,N);
Y = 1/N*(X-repmat(mean(X),[N,1]))'*(X-repmat(mean(X),[N,1]));
obj = clime(Y, 1e-7);
omega = getinstance(obj,1e-7);
fprintf('Error measured in L2 norm is %12.8f\n', norm(omega-inv(Y)));
fprintf('Error measured in infinity norm is %12.8f\n', norm(omega-inv(Y),Inf));

%%
%Example of clime estimator. Same setting with previous example, but this time 
%n < p.
clear all
m = 20;
N = 20;
S=diag(ones(2*m+1,1)) + 0.5*diag(ones(2*m,1),1) + 0.5*diag(ones(2*m,1),-1);
X = mvnrnd(zeros(1,2*m+1), S,N);
Y = 1/N*(X-repmat(mean(X),[N,1]))'*(X-repmat(mean(X),[N,1]));
obj = clime(Y, 1e-7,1e-4);
omega = getinstance(obj,1e-5);
fprintf('Error measured in L2 norm is %12.8f\n', norm(omega*Y-eye(size(Y))));
fprintf('Error measured in infinity norm is %12.8f\n', max((max(omega*Y-eye(size(Y))))));

%%
%Exapmle of tunning parameter, once we have the solution path object, it's
%very convenient to tune the parameters. 
X = rand(50,100);
beta = sprand(100,1,0.1)*10;
Y = X*beta + rand(50,1)*0.1;
obj = dantzig(X,Y);
soln = getinstance(obj,0.1);
fprintf('Choosing lambda = 1e-1, error in infinity norm is %12.8f\n',...
norm(soln-beta, Inf));
soln = getinstance(obj, 1e-2);
fprintf('Choosing lambda = 1e-2, error in infinity norm is %12.8f\n',...
norm(soln-beta,Inf));
soln = getinstance(obj, 1e-3);
fprintf('Choosing lambda = 1e-3, error in infinity norm is %12.8f\n',...
norm(soln-beta,Inf));
soln = getinstance(obj, 1e-4);
fprintf('Choosing lambda = 1e-4, error in infinity norm is %12.8f\n',...
norm(soln-beta,Inf));








