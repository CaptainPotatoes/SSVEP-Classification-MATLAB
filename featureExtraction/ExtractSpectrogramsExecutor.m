%% Extract Spectrograms % SSVEP
clear;clc;close all;
% [DATA,filename] = csvread('Subject1_SingleChannel_10Hz_to_16Hz.csv');
% [DATA,filename] = csvread('EEG_SSVEPData_2017.05.31_14.55.24.csv');
[DATA, filename] = csvread('Subject1_Trial1.1.csv');
Fs = 250;
X_1 = DATA(:,1);
X_2 = DATA(:,2);
h = 1/250;
t=0:h:(size(DATA,1)/250)-h;
range = 250:250:1000;
start = 1;
wStart = start:250:(length(X_1)-max(range));
PLOTDATA = 1==0;
%% Filt/Extract Spectrograms %{
for i = 1:length(wStart)
    start = (i-1)*115+1; 
    fin = start + 114;%S(start:fin,:)
    [T,F,STFT1{i}] = extractSpectrograms(X_1(wStart(i):wStart(i)+999),PLOTDATA);
    [~,~,STFT2{i}] = extractSpectrograms(X_2(wStart(i):wStart(i)+999),PLOTDATA);
    CLASS{i} = unique(DATA(wStart(i):wStart(i)+999,3));
%     a = input('Continue? \n');
end
clearvars -except filename T F STFT1 STFT2 CLASS
save([filename(1:end-4) '_spect.mat'],'-v7','T','F','STFT1','STFT2','CLASS')
%}