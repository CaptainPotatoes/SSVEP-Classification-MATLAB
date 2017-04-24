function [ y ] = wk1mod( x )
%UNTITLED4 mod of wkeep1
len = 250;
y = x;
sx = length(x);
ok = (len >= 0) && (len < sx);
if ~ok , return; end
d = (sx-len)/2;
first = 1+floor(d); last = sx-ceil(d);
y = y(first:last);

end

