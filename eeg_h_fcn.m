function [ecg_h] = eeg_h_fcn(ecg, fs)
%ECG_H_FCN Summary of this function goes here
%   Detailed explanation goes here
if ~isvector(ecg)
  error('ecg must be a row or column vector');
end
ecg = ecg(:); %vectorize

f2s = fs;
%% bandpass filter for Noise cancelation of other sampling frequencies(Filtering)
f1=13; %cuttoff low frequency to get rid of baseline wander
f2=90; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/f2s; % cut off based on fs
N = 3; % order of 3 less processing
[a,b] = butter(N,Wn); %bandpass filtering
ecg_h_a = filtfilt(a,b,ecg);
ecg_h_a = ecg_h_a/max(abs(ecg_h_a));
%% Notch at 60Hz
cc = [59 61];
c0 = cc/(f2s/2);
[b1, a1] = butter(3,c0, 'stop');
ecg_h_60 = filtfilt(b1,a1,ecg_h_a);
cc2 = [49 51];
c02 = cc2/(f2s/2);
[b2, a2] = butter(3,c02, 'stop');
ecg_h = filtfilt(b2,a2,ecg_h_60);
% d = designfilt('bandstopiir','FilterOrder', 2, ...
%                'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
%                'DesignMethod','butter','SampleRate',fs);

% d3 = designfilt('bandstopiir','FilterOrder', 2, ...
%                'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
%                'DesignMethod','butter','SampleRate',fs);


%% Notch filter at 50Hz;
% d3 = designfilt('bandstopiir','FilterOrder',2, ...
%                'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
%                'DesignMethod','butter','SampleRate',fs);
% fvtool(d,'Fs',fs)
% ecg_h_60 = filtfilt(d, ecg_h_a);
% ecg_h_c = filtfilt(d3, ecg_h_b);
% ecg_h_d = filtfilt(d4, ecg_h_c);
% ecg_h = filtfilt(d3, ecg_h_60);
end