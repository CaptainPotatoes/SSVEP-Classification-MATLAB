function [ecg_h] = eeg_h_fcn(ecg, fs)
%ECG_H_FCN Summary of this function goes here
%   Detailed explanation goes here
if ~isvector(ecg)
  error('ecg must be a row or column vector');
end
ecg = ecg(:); %vectorize
f2s = fs;
%% bandpass filter for Noise cancelation of other sampling frequencies(Filtering)
f1=3; %cuttoff low frequency to get rid of baseline wander
f2=50; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/f2s; % cut off based on fs
N = 3; % order of 3 less processing
[a,b] = butter(N,Wn); %bandpass filtering
ecg_h_a = filtfilt(a,b,ecg);
% ecg_h_a = ecg_h_a/max(abs(ecg_h_a));
%% Notch at 60Hz
% cc = [59 61];
% c0 = cc/(f2s/2);
% [b1, a1] = butter(3,c0, 'stop');
% ecg_h = filtfilt(b1,a1,ecg_h_a);

cc2 = [49 51];
c1 = cc2/(f2s/2);
[b2, a2] = butter(3,c1, 'stop');
ecg_h = filtfilt(b2,a2,ecg_h_a);
end