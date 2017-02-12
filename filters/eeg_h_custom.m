function [Y] = eeg_h_custom(X,fs,f,N)
%eeg_h_custom
if ~isvector(X)
  error('ecg must be a row or column vector');
end
X = X(:); %vectorize
f2s = fs;
%% bandpass filter for Noise cancelation of other sampling frequencies(Filtering)
% f1=4; %cuttoff low frequency to get rid of baseline wander
% f2=50; %cuttoff frequency to discard high frequency noise
Wn=[f(1) f(2)]*2/f2s; % cut off based on fs
% N = 7; % order of 3 less processing
[a,b] = butter(N,Wn); %bandpass filtering
Y1 = filtfilt(a,b,X);
% ecg_h_a = ecg_h_a/max(abs(ecg_h_a));

cc2 = [49 51];
c1 = cc2/(f2s/2);
[b2, a2] = butter(3,c1, 'stop');
Y = filtfilt(b2,a2,Y1);
end