function [ f, C ] = get_fft_data( X, Fs )
%get_fft_data:
% X is filtered data
L = size(X,1);
A = fft(X);
B = abs(A/L);
C = B(1:(L/2+1));
C(2:end-1) = 2*C(2:end-1);
f = Fs*(0:(L/2))/L;

end

 