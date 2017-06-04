% Experiments
%% Generating Reference Signal:
clear;close all;clc
ic = 0.1; C1 = 9.7:ic:10.3; C2 = 12.1:ic:12.7; C3 = 14.8:ic:15.4; C4 = 16.2:ic:16.8;
f_new = [C1,C2,C3,C4]; Fs = 250;F = [4 37]; winLim = [5 37];
amplitude = 0.5E-4;
[DATA, filename] = csvread('Subject1_Trial1.1.csv');
numch = 2; datach = DATA(:,1:numch);
filtch = zeros(size(datach,1),numch);
for i = 1:numch %    
	filtch(:,i) = customFilt(datach(:,i),Fs,F,3); %figure(1); hold on; plot(filtch(:,i));
end
for i = 1:length(f_new)
    [sigs(i,:)] = testSignal(f_new(i),5000,amplitude,250);
end
sumSigs = sum(sigs);
[S1,wfreqs] = welch_psd(sumSigs, Fs, hannWin(2048));
plot(wfreqs, S1),xlim(winLim);xlabel('Frequency (Hz)'),ylabel('Power Density (W/Hz)'),title([ filename, ' - ',  'Power Spectral Density Estimate']);
