function [ Y ] = scaleMin( X )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
X1 = sort (X);
X2 = X1(X1~=0);
if ~isempty(X2)
    Y = X/X2(1);
else
    Y = X;
end
end

