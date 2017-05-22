%% 
clear;close all;clc

datach = csvread('Matt_1ch_10_to_16_3.csv');
hannWin = hann(2048); wlen = 1024; h=64; nfft = 4096;
filtch = zeros(size(datach,1),1);
i=1; Fs = 250;

filtch(:,i) = customFilt(datach(:,i),Fs,[8 20],3); %figure(1); hold on; plot(filtch(:,i));

for i = 1:1000
%     [S1, f1, t1] = stft( filtch(:,1), wlen, h, nfft, Fs );
    [S2, f2, t2] = stft2( filtch(:,1), wlen, h, nfft, Fs );
end