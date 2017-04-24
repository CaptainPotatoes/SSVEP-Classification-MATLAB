function [ y ] = wc1mod( x,f )
%UNTITLED5 mod of wconv1
y = conv2(x(:)',f(:)','full');
if size(x,1)>1
    y = y';
    
end

