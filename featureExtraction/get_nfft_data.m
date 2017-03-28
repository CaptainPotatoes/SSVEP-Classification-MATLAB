function [ f, C ] = get_nfft_data( X, Fs, L )
%get_fft_data:
% X is filtered data
% L = size(X,1);
% L = number of FFT points
L = 2^nextpow2(L);
A = fft(X,L);
B = abs(A/L);
%% TEMP : ONLY FOR EVEN SIZED ARRAYS:
% C = zeros(1,(L/2+1));
C = B(1:(L/2+1));
C(2:end-1) = 2*C(2:end-1);
%{
if mod(L,2) == 0
    C = B(1:(L/2+1));
    C(2:end-1) = 2*C(2:end-1);
elseif mod(L,2) == 1
    C = B(1:(L/2+0.5));
    C(2:end-1) = 2*C(2:end-1);
end
%}
f = Fs*(0:(L/2))/L;

end
