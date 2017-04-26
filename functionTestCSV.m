clear;clc;close all;
% LOAD TRAINING DATA: (tX, tY);
datach = csvread('FadiTest1.csv');
removeStart = 0;
removeEnd   = 0;
ch1 = datach(:,1);
ch2 = datach(:,2);
ch3 = datach(:,3);
Fs = 250;
% samplesX = ch1(end-249:end);
% waveletCoef1 = cwt(samplesX,1:20,'haar');
% waveletCoef2 = cwt_haar(samplesX);
%%-Plot Analysis: 
%{
winLim = [6 24];
filtch = zeros(size(datach,1),size(datach,2));
hannWin = hann(4096); 
wlen = 1024; h=64; nfft = 4096;
K = sum(hamming(wlen, 'periodic'))/wlen;
figure(1);plot(datach);
figure(2);hold on;%FFT
for i = 1:4
    filtch(:,i) = eegcfilt(datach(:,i)); plot(filtch(:,i));
%     [f, P1] = get_fft_data(filtch(:,i),Fs);
%     plot(f,P1),xlim(winLim);
end
figure(3);hold on;%PSD
for i = 1:1
    [S1,wfreqs] = welch_psd(filtch(:,i), 250, hannWin); 
    plot(wfreqs, S1),xlim(winLim);
end   
fH = figure(4);hold on; set(fH, 'Position', [-2560, 0, 1600, 900]);%Spect
for i = 1:1
    subplot(2,2,i);
    [S1, f1, t1] = stft( filtch(:,i), wlen, h, nfft, Fs ); S2 = 20*log10(abs(S1(f1<winLim(2) & f1>winLim(1),:))/wlen/K + 1e-6); 
    imagesc(t1,f1(f1<winLim(2) & f1>winLim(1)),S2);set(gca,'YDir','normal');xlabel('Time, s');ylabel('Frequency, Hz');colormap(jet)
    cb = colorbar;
    ylabel(cb, 'Power (db)')
    title(['Ch' num2str(i)]);
end
%}
%% Classification EOG;
winSize = 250; %1s:
c=1;
for i = 1:250:length(ch1)-250
    ch1_p = ch1(i:i+249);
    ch2_p = ch2(i:i+249);
    ch3_p = ch3(i:i+249);
    Y(c) = eog_knn(ch1_p,ch2_p,ch3_p);
    c=c+1;
end


%% Classification:
%{
range = 500:60:2500;
Window = cell(size(range,2),4);
Y = cell(size(range,2),1);
cont = [];
EOGONLY = false;
PLOTDATA = isempty(cont);
OUT = zeros(1,size(range,2));
History = zeros(size(range,2),4);
for i = 1:size(range,2)
    start = 1;
    winEnd = start + (range(i)-1);
    fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    fprintf('Seconds Elapsed = [%1.2f]\r\n',winEnd/250);
    Window{i,1} = ch1( start : winEnd ); % set values:
    Window{i,2} = ch2( start : winEnd );
    Window{i,3} = ch3( start : winEnd );
    Window{i,4} = ch4( start : winEnd );
    [Y{i},F{i}] = fHC(Window{i,1}, Window{i,2}, Window{i,3}, ...
        Window{i,4}, Fs, EOGONLY, PLOTDATA);
    [History(i,:), OUT(i)] = featureAnalysis(F{i},winEnd);
    meanH = mean(History(1:i,:))
    if OUT(i)~=0
        countH(i) = countOccurrences(OUT(:,1:i), OUT(i));
    else
        countH(i) = 0;
    end
    if (max(meanH)>7) && countH(i)>=5
        OUTPROPER(i) = OUT(i);
    else
        OUTPROPER(i) = 0;
    end
    if isempty(cont)
        commandwindow;
        cont = input('Continue? \n');
    end
end
%}
% % EOF
