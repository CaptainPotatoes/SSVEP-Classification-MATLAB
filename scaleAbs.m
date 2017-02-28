function [ Y ] = scaleAbs( X )
%scaleAbs :: scales to the maximum of the absolute value of the signal:

X = X(:);
xmax = max(abs(X));
Y = X/xmax;

end

