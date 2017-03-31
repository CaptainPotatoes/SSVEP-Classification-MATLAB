function [ X ] = scaleMax( X )
F = find(X);
X(F)=X(F)-min(X);
end

