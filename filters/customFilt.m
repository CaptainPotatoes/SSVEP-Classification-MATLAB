function [ Y ] = customFilt( X,Fs,f,N )
%FILT_CUSTOM Allows for quick customization of bandpass filter parameters
if ~isvector(X)
  error('must be a row or column vector');
end
X = X(:); %vectorize
Wn=[f(1) f(2)]*2/Fs; % cut off based on Fs
[a,b] = butter(N,Wn);
Y = filtfilt(a,b,X);
end

