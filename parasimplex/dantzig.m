% dantzig.m
% dantzig solves the constrained L1 minimization problem:
%   min ||b||_1
%   s.t. 
%     ||X'*(Y-Xb)||_inf <= lambda
%
% This problem is recast as a linear program:
%   min b'+b
%   s.t.
%     X'*(Y-X(b'-b)) <= lambda
%    -X'*(Y-X(b'-b)) <= lambda
% and is solved using parametric simplex method. The parametric simplex
% method is a homotopy method, it starts with lambda=inf, then decreases
% lambda along the way to zero. Solutions corresponding to all lambda values
% are computed at once.
%
% Usage: soln = dantzig(X, Y, lambda); 
%
% Input parameters:
% 
% X  -- n*p matrix, the same as defined in the minimization problem.
%
% Y  -- n*1 vector, the same as defined in the minimization problem.
%
% lambda -- scalar(optional), the same as defined in the minimizaiton
%           problem. By default, lambda is set to be 1e-6.
%
% Output parameters:
%
% soln -- a struct which contains the solution path. Can be used by function
%         getinstance to retrieve the solution to a specific lambda value.
%         See getinstance.m for more details.
%
% Written by: Cheng Tai, Princeton.
% Parametric simplex solver written by Robert Vanderbei, Princeton.
% Email: chengt@princeton.edu
% Date: January, 2013.

function [ soln ] = dantzig(X, Y, lambda)
if nargin < 2
    error('Too few input arguments.');
elseif nargin > 3
    error('Too many input arguments.');
elseif nargin == 2
    lambda = 1e-6;
end

if size(X,1) ~= size(Y,1) || size(Y,2) ~= 1
    error('Input format not correct.');
end

%set up of the problem
disp('Setting up problems...');
p = size(X,2);
T = X'*X;
A = [T -T; -T T];
b1 = X'*Y;
b = [b1; -b1];
c = -ones(2*p,1);
disp('Calling parametric simplex solver...');
%plug in the solver
[R M] = mex_dantzig(A,b,c,lambda);
soln.type = 'dantzig';
soln.R = R;
soln.M = M;

disp('Computation completed.');

end
