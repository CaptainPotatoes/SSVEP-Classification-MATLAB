function [eog_h] = eog_h_fcn(eog, fs)
%ECG_H_FCN Summary of this function goes here
%   Detailed explanation goes here
if ~isvector(eog)
  error('ecg must be a row or column vector');
end
eog = eog(:); %vectorize
f2s = fs;
%% bandpass filter for Noise cancelation of other sampling frequencies(Filtering)
fz=10; %cuttoff frequency to discard high frequency noise
Wn=fz*2/f2s; % cut off based on fs
N = 3; % order of 3 less processing
[a,b] = butter(N,Wn,'low'); %bandpass filtering
eog_h = filtfilt(a,b,eog);
%% Temp
fh=2; %cuttoff frequency to discard high frequency noise
Wn=fh*2/f2s; % cut off based on fs
N = 3; % order of 3 less processing
[a,b] = butter(N,Wn,'high'); %bandpass filtering
eog_h = filtfilt(a,b,eog_h);

end