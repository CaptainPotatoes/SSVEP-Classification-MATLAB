%% Filter Data:
clear;clc;close all;
[DATA,filename] = csvread('Subject1_Trial1.1.csv');
% Fs = 250;
freqs = [250,500,1000,2000,4000];
windowSize = 4*freqs;
datach = DATA(:,1:2);
% classLabels = DATA(:,3);
F = [4 35];winLim = [8 35];
filtch = cell(length(freqs),1);
for i = 1:length(freqs)
    for wins = 1:length(windowSize)
        filtch{i} = ecg_var_filter(datach(1:windowSize(i),1), freqs(i), windowSize(i));
    end
%    filtch(:,i) = customFilt(datach(:,i),Fs,F,3);
%     filtch(:,i) = ssvep_filter_f32(datach(1:windowSize,i));
end
%% Variable Window Sizes:
fH = figure(4);
set(fH, 'Position', [0, 0, 1200, 1400]);%Spect
plot(filtch{1})
for i = 1:length(freqs)
    windowSize = [1:8]*freqs(i);
    for wins = 1:length(windowSize)
        filtch_varsize{i,wins} = ecg_var_filter(datach(1:windowSize(wins),1), freqs(i), windowSize(wins));
    end
%    filtch(:,i) = customFilt(datach(:,i),Fs,F,3);
%     filtch(:,i) = ssvep_filter_f32(datach(1:windowSize,i));
end
%{
wlen = 1024; h=64; nfft = 4096;
K = sum(hamming(wlen, 'periodic'))/wlen;
for i = 1:2
    subplot(2,1,i)
    [S1, f1, t1] = stft2( filtch(:,i), wlen, h, nfft, Fs ); S2 = 20*log10(abs(S1(f1<winLim(2) & f1>winLim(1),:))/wlen/K + 1e-6); 
    imagesc(t1,f1(f1<winLim(2) & f1>winLim(1)),S2),xlim([min(t1) max(t1)]),ylim(winLim);
    set(gca,'YDir','normal');xlabel('Time, s');ylabel('Frequency, Hz');colormap(jet)
    cb = colorbar;ylabel(cb, 'Power (db)')
    title(['Ch' num2str(i)]);
end
%}
% DATAWRITE = [filtch, classLabels];
% csvwrite([filename(1:end-4) '_filt.csv'],DATAWRITE);