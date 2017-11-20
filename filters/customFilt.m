function [ Y, b, a ] = customFilt( X,Fs,f,N )
%FILT_CUSTOM Allows for quick customization of bandpass filter parameters
if ~isvector(X)
  error('must be a row or column vector');
end
% coder.varsize('a','b');
X = X(:); %vectorize
Wn=[f(1) f(2)]*2/Fs; % cut off based on Fs
[b,a] = butter(N,Wn);
Y = filtfilt(b,a,X);
end

