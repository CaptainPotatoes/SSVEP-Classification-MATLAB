function [ Y ] = scaleAbs( X )
%scaleAbs :: scales to the maximum of the absolute value of the signal ...
% (For a ratio of 1).

X = X(:);
xmax = max(abs(X));
Y = X/xmax;

end

