function [ w ] = hammPeriodic( x )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x=x+1;
if ~rem(x,2)
    % Even length window
    half = x/2;
    w = calc_window(half,x);
    w = [w; w(end:-1:1)];
else
    % Odd length window
    half = (x+1)/2;
    w = calc_window(half,x);
    w = [w; w(end-1:-1:1)];
end
function w = calc_window(m,x)
    c = (0:m-1)'/(x-1);
    w = 0.54 - 0.46*cos(2*pi*c);
end
end

