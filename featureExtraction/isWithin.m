function [ B ] = isWithin( X, range )
%isWithin Summary of this function goes here
%   Checks if variable 'X' is within frange
% X is a scalar
% range is a 2x1 or 1x2 vector, with min and max values respectively
% Check Range
if(size(range,1)*size(range,2) ~= 2)
    error('ERROR: Incorrect Range!');
end
% Vectorize range:
range = range(:); 
B = X>range(1,1) && X<range(2,1);

end

