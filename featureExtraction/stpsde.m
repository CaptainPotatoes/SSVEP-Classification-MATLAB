function [ S,F,T ] = stpsde( x, wlen, h, nfft, fs )
%% As of 5/26 - does not work yet. 
%STPSDE - Experimental Short-term Power Spectral Density Estimation
% x - signal in the time domain
% wlen - length of the hamming window
% h - hop size
% nfft - number of FFT points
% fs - sampling frequency, Hz
% F - frequency vector, Hz
% T - time vector, S
% S - STFT matrix (only unique points, time across columns, freq across rows)
x=x(:);
xlen = length(x);

win = hammPeriodic(wlen);
rown = ceil((1+nfft)/2);            % calculate the total number of rows
coln = 1+fix((xlen-wlen)/h);        % calculate the total number of columns
% S = complex(zeros(rown, coln));           % form the stft mat
hannW = hannWin(256);
indx = 0:h:(xlen-wlen);
for c = 1:size(indx,2)
    xw = x(indx(c)+1:indx(c)+wlen).*win;
    X = welch_psd2(xw,fs,hannW);
    S(:, c) = X(1:rown-1);
end
T = (wlen/2:h:wlen/2+(coln-1)*h)/fs;
F = (0:rown-1)*fs/nfft;
K = sum(hammPeriodic(wlen))/wlen;
figure(7);imagesc(T,F,(20*log10(abs(S))/wlen/K+1E-6));
end

