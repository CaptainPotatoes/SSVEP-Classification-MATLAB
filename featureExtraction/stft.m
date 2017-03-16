%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Short-Time Fourier Transform            %
%               with MATLAB Implementation             %
%                                                      %
% Author: M.Sc. Eng. Hristo Zhivomirov       12/21/13  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s, f, t] = stft(x, wlen, h, nfft, fs)

% function: [s, f, t] = stft(x, wlen, h, nfft, fs)
% x - signal in the time domain
% wlen - length of the hamming window
% h - hop size
% nfft - number of FFT points
% fs - sampling frequency, Hz
% f - frequency vector, Hz
% t - time vector, s
% s - STFT matrix (only unique points, time across columns, freq across rows)

% represent x as column-vector if it is notx = x(:);
if size(x, 2) > 1
    x = x';
end

% length of the signal
xlen = length(x);
% wlen = 2^nextpow2(fs);
% form a periodic hamming window (!REPLACED EXTERNALLY)
% win = hamming(256, 'periodic');
win = hammPeriodic(wlen);
%WINDOW LENGTH:


% form the stft matrix
rown = ceil((1+nfft)/2);            % calculate the total number of rows
coln = 1+fix((xlen-wlen)/h);        % calculate the total number of columns
s = complex(zeros(rown, coln));           % form the stft matrix
% for r = 1:rown
%     for c = 1:coln
%         s(r,c) = complex(0,0);
%     end
% end
% initialize the indexes
indx = 0;
col = 1;

% perform STFT
while indx + wlen <= xlen
    % windowing
    xw = x(indx+1:indx+wlen).*win;
    
    % FFT
    X = fft(xw, nfft);
    
    % update the stft matrix
%     s(:, col) = X(1:rown); %%%% OLD
    s(:, col) = X(1:rown);
    
    % update the indexes
    indx = indx + h;
    col = col + 1;
end

% calculate the time and frequency vectors
t = (wlen/2:h:wlen/2+(coln-1)*h)/fs;
f = (0:rown-1)*fs/nfft;

end