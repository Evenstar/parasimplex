% getinstance.m
% This is a utility function for retrieving a solution from a solution path
% struct. As the solution path for dantzig and clime is piecewise linear, this
% function get the solution corresponding to a specific lambda value by linear
% interpolation.
%
% Usage: soln = getinstance(obj, lambda);
%
% Input parameters:
%
% obj  -- solution path struct, returned by calling function dantzig or
%         clime.
%
% lambda -- scalar, corresponding to which value the solution is
%           interpolated. Note this value should be not exceed the range of
%           lambda values of obj. For example, if obj is computed until
%           1e-6, then this lambda can not be smaller than 1e-6.
%
% Output parameters:
%
% soln -- vector or matrix. For dantzig, soln is a p*1 vector. For clime,
%         soln is p*p matrix, note for clime, the solution of the linear
%         program is usually not symmetric, so soln is symmetrized by
%         setting soln(i,j)=soln(j,i)=the one with smaller absolute value
%         between soln(i,j) and soln(j,i).
%
% Written by: Cheng Tai, Princeton
% Email: chengt@princeton.edu
% Date: January, 2013
function [ soln ] = getinstance(obj, lambda)
if ~isstruct(obj)
    error('First input argument is not a struct.');
end
if strcmp(obj.type,'dantzig')
    soln = interp1(obj.R, obj.M, lambda);
    soln = soln';
end
if strcmp(obj.type, 'clime')
    soln = zeros(size(obj.cR{1},2));
    for i = 1 : size(obj.cR,1)
        beta_icol = interp1(obj.cM{i}, obj.cR{i}, lambda);
        soln(:,i) = beta_icol;
    end
    %Symmetrization
    soln = (abs(soln) <= abs(soln')).*soln + (abs(soln) > abs(soln')).*soln';
end
end
