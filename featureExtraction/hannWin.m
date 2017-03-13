function [ w ] = hannWin( x )
%Calculate the generalized cosine window samples
% x is the length of the window

w = .5*(1 - cos(2*pi*(0:x-1)'/(x-1)));

end

