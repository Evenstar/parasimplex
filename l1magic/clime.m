% clime.m
%
% clime is used to estimate the precision matrix(inverse of covariance
% matrix). It solves the constrained L1 minimization problem:
%   min ||S||_1
%   s.t.
%     ||(X+rI)*S - I||_inf <= lambda
%
% This problem is recast as a linear program and is solved in a column by
% column fashion.
%
% For each i in 1..p, p is dimension of X, solve
%   min ||b_i||_1
%   s.t.
%     ||(X+rI)*b_i -e_i||_inf <= lambda
% where b_i is the i-th column of S, and e_i is the i-th unit vector with
% e(i)=1.
%  
% Then this family of problems is solved by parametric simplex method.
%
% Usage: obj = clime(X, lambda, r);
%
% Input parameters:
% 
% X -- p*p matrix, X is sample covariance matrix.
%
% lambda -- scalar(optional), the same as defined in the problem. By
%           default, lambda is set to be 1e-6.
%
% r -- scalar, the same as defined in the problem. This is originally
%      introduced for numerical stability.
%
% Output parameters:
%
% obj -- a struct contains the solution path. Can be used by function
%        getinstance to retrieve the solution to a specific value of
%        lambda. See getinstance.m for more details.
%
% Written by: Cheng Tai, Princeton.
% Parametric simplex solver written by Robert Vanderbei, Princeton.
% Email: chengt@princeton.edu
% Date: January, 2013.
%
function [ obj ] = clime(X, lambda, r)
if nargin < 1
    error('Too few input arguments.');
elseif nargin > 3
    error('Too many input arguments.');
elseif nargin == 1
    lambda = 1e-6;
    r = 0;
elseif nargin == 2
    r = 0;
else 
   X = X + r*eye(size(X));
end

if size(X,1) ~= size(X,2)
    error('Input format not correct.');
end

A = [X -X; -X X];
disp('Calling parametric simplex solver...');
[cR cM] = mex_clime(A, lambda);

obj.type = 'clime';
obj.cR = cR;
obj.cM = cM;

disp('Computation completed.');
end

