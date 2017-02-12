function [ f, C ] = get_fft_data( X, Fs )
%get_fft_data:
% X is filtered data
L = size(X,1);
A = fft(X);
B = abs(A/L);
if mod(L,2) == 0
    C = B(1:(L/2+1));
    C(2:end-1) = 2*C(2:end-1);
elseif mod(L,2) == 1
    C = B(1:(L/2+0.5));
    C(2:end-1) = 2*C(2:end-1);
end
f = Fs*(0:(L/2))/L;

end

 