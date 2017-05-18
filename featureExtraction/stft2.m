%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Short-Time Fourier Transform            %
%               with MATLAB Implementation             %
%                                                      %
% Author: M.Sc. Eng. Hristo Zhivomirov       12/21/13  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S, F, T] = stft2(x, wlen, h, nfft, fs)

% function: [S, F, T] = stft2(x, wlen, h, nfft, fs)
% x - signal in the time domain
% wlen - length of the hamming window
% h - hop size
% nfft - number of FFT points
% fs - sampling frequency, Hz
% F - frequency vector, Hz
% T - time vector, S
% S - STFT matrix (only unique points, time across columns, freq across rows)

% represent x as column-vector if it is not x = x(:);
% if size(x, 2) > 1
%     x = x';
% end
x=x(:);
% length of the signal
xlen = length(x);
win = hammPeriodic(wlen);
% form the stft matrix
rown = ceil((1+nfft)/2);            % calculate the total number of rows
coln = 1+fix((xlen-wlen)/h);        % calculate the total number of columns
S = complex(zeros(rown, coln));           % form the stft matrix
% initialize the indexes
% perform STFT
%{
while indx + wlen <= xlen
    % windowing
    xw = x(indx+1:indx+wlen).*win;
    
    % FFT
    X = fft(xw, nfft);
    
    % update the stft matrix
    S(:, col) = X(1:rown);
    
    % update the indexes
    indx = indx + h;
    col = col + 1;
end
%}
%Correct way to perform STFT:
indx = 0:h:(xlen-wlen);
for c = 1:size(indx,2)
    xw = x(indx(c)+1:indx(c)+wlen).*win;
    X = fft(xw, nfft);
    S(:, c) = X(1:rown);
end
% calculate the time and frequency vectors
T = (wlen/2:h:wlen/2+(coln-1)*h)/fs;
F = (0:rown-1)*fs/nfft;
% Select Relevant Part of Signal:

end