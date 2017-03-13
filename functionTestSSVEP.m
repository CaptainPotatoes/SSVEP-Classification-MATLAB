    %Clear 
clear;close all;clc;
    %Import Data:
ChannelNames = {['Fp1' 'Fp2' 'Fpz' 'REye']};
% load('mssvep_16.6_3.mat');
% load('mssvep_15_1.mat');
load('mssvep_10_2.mat');
% load('mssvep_12.5_1.mat')
% load('mssvep_t1_baseline');
% load('mssvep_t2_16_1.mat');
%-SHOW SPECT:
showSpect = 0;
%--- LOAD CLASS ---%
CLASS = '10';
VERSION = 'v1';
%--- START ANALYSIS (PART1) ---%
remove = 0; % Remove final second of data.
removeFromStart = 0;

Fs = SamplingRate;
%Import as variables and scale all to one:
ch1 = Trial{1}(1+removeFromStart:end-remove,1);
ch2 = Trial{2}(1+removeFromStart:end-remove,1);
ch3 = Trial{3}(1+removeFromStart:end-remove,1);
ln = min([length(ch1) length(ch2) length(ch3)]);

ch1 = ch1(1:ln);
ch2 = ch2(1:ln);
ch3 = ch3(1:ln);

if size(Trial,2) > 3
    ch4 = Trial{4}(1:end-remove,1);
end

seconds = length(ch1)/Fs
flim   = [8.0 18];
winLim = [9 17.3];
%
N = 5;
    %Filter & Scale everything to '1'
ch1_f = scaleAbs(customFilt(ch1, Fs, flim, N));
ch2_f = scaleAbs(customFilt(ch2, Fs, flim, N));
ch3_f = scaleAbs(customFilt(ch3, Fs, flim, N));
if size(Trial,2) > 3
    ch4_f = scaleAbs(customFilt(ch4, Fs, flim, N));
end
[f,  P1]  = get_fft_data(ch1_f, Fs);
[~, P2] = get_fft_data(ch2_f, Fs);
[~, P3] = get_fft_data(ch3_f, Fs);
if size(Trial,2) > 3
    [f4, P4] = get_fft_data(ch4_f, Fs);
end
    fH = figure(1); 
    set(fH, 'Position', [100, 100, 1440, 820]);
    subplot(2,2,1);
hold on;
plot(f,  P1,'color','m'),xlim([1 35]);
plot(f, P2,'color','c'),xlim([1 35]);
plot(f, P3,'color','r'),xlim([1 35]);
if size(Trial,2) > 3
    plot(f4, P4,'color','b'),xlim([1 35]);
end
hold off;
title('FFT(Ch1-4)');
ylabel('|P1(f)|');
xlabel('f (Hz)');
    subplot(2,2,2)
    plot(f,(P1+P2+P3)),xlim([1 35]);
subplot(2,2,3);
% wind = [1024 512 256 128];
hannWin = hann(1024);
[S1,wfreqs] = welch_psd(ch1_f, 250, hannWin); 
[S2,~    ] = welch_psd(ch2_f, 250, hannWin);
[S3,~    ] = welch_psd(ch3_f, 250, hannWin);
hold on;
plot(wfreqs, S1),xlim([1 35]);
plot(wfreqs, S2),xlim([1 35]);
plot(wfreqs, S3),xlim([1 35]);
hold off;

subplot(2,2,4)
plot(wfreqs, (S1+S2+S3)),xlim([1 35]);
% Spectrograms:
if showSpect == 1
    f2 = figure(2);
    set(f2, 'Position', [100, 100, 1200, 675]);
        subplot(2,2,1)
    [~, Fspect, T, P] = spectrogram(ch1_f, 5*Fs,4*Fs,10*Fs,Fs);
    imagesc(T, Fspect(Fspect<winLim(2) & Fspect>winLim(1)), 10*log10(P(Fspect<winLim(2) & Fspect>winLim(1),:)));
    set(gca,'YDir','normal')
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    cb = colorbar;
    ylabel(cb, 'Power (db)')
    colormap(jet)
    title('Channel 1', 'FontSize', 14)
        subplot(2,2,2)
    [~, Fspect, T, P2] = spectrogram(ch2_f, 5*Fs,4*Fs,10*Fs,Fs);
    imagesc(T, Fspect(Fspect<winLim(2) & Fspect>winLim(1)), 10*log10(P2(Fspect<winLim(2) & Fspect>winLim(1),:)));
    set(gca,'YDir','normal')
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    cb = colorbar;
    ylabel(cb, 'Power (db)')
    colormap(jet)
    title('Channel 2', 'FontSize', 14)
        subplot(2,2,3)
    [~, Fspect, T, P3] = spectrogram(ch3_f, 5*Fs,4*Fs,10*Fs,Fs);
    imagesc(T, Fspect(Fspect<winLim(2) & Fspect>winLim(1)), 10*log10(P3(Fspect<winLim(2) & Fspect>winLim(1),:)));
    set(gca,'YDir','normal')
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    cb = colorbar;
    ylabel(cb, 'Power (db)')
    colormap(jet)
    title('Channel 3', 'FontSize', 14)
        subplot(2,2,4)
        P1_3combined =  10*log10(P(Fspect<winLim(2) & Fspect>winLim(1),:))+10*log10(P2(Fspect<winLim(2) & Fspect>winLim(1),:))...
            +10*log10(P3(Fspect<winLim(2) & Fspect>winLim(1),:));
    imagesc(T, Fspect(Fspect<winLim(2) & Fspect>winLim(1)),P1_3combined);
    set(gca,'YDir','normal')
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    cb = colorbar;
    ylabel(cb, 'Power (db)')
    colormap(jet)
    title('Channel Sum(1-3)', 'FontSize', 14)
   
    f3 = figure(3);
    set(f3, 'Position', [100, 100, 1200, 675]);
    % wlen = 2^nextpow2(Fs);
    wlen = 4*Fs;
    h=256;
    nfft = 4096;
    K = sum(hamming(wlen, 'periodic'))/wlen;
        subplot(2,2,1)
    [s1, f1, t1] = stftOrig( ch1_f, wlen, h, nfft, Fs );
    s1_1 = 20*log10(abs(s1(f1<winLim(2) & f1>winLim(1),:))/wlen/K + 1e-6); 
    imagesc(t1,f1(f1<winLim(2) & f1>winLim(1)),s1_1);
    set(gca,'YDir','normal')
    xlabel('Time, s')
    ylabel('Frequency, Hz')
    title('Amplitude spectrogram of Ch1')
    handl = colorbar;
    colormap(bone)
    ylabel(handl, 'Magnitude, dB')
        subplot(2,2,2)
    [s2, ~, ~] = stftOrig( ch2_f, wlen, h, nfft, Fs );
    s2_1 = 20*log10(abs(s2(f1<winLim(2) & f1>winLim(1),:))/wlen/K + 1e-6); 
    imagesc(t1,f1(f1<winLim(2) & f1>winLim(1)),s2_1);
    set(gca,'YDir','normal')
    xlabel('Time, s')
    ylabel('Frequency, Hz')
    title('Amplitude spectrogram of Ch2')
    handl = colorbar;
    colormap(bone)
    ylabel(handl, 'Magnitude, dB')
        subplot(2,2,3)
    [s3, ~, ~] = stftOrig( ch3_f, wlen, h, nfft, Fs );
    s3_1 = 20*log10(abs(s3(f1<winLim(2) & f1>winLim(1),:))/wlen/K + 1e-6); 
    imagesc(t1,f1(f1<winLim(2) & f1>winLim(1)),s3_1);
    set(gca,'YDir','normal')
    xlabel('Time, s')
    ylabel('Frequency, Hz')
    title('Amplitude spectrogram of Ch3')
    handl = colorbar;
    colormap(bone)
    ylabel(handl, 'Magnitude, dB')
        subplot(2,2,4)
    s4_1 = 20*log10(abs(s1(f1<winLim(2) & f1>winLim(1),:))/wlen/K + 1e-6)+...
        20*log10(abs(s2(f1<winLim(2) & f1>winLim(1),:))/wlen/K + 1e-6)+...
        20*log10(abs(s3(f1<winLim(2) & f1>winLim(1),:))/wlen/K + 1e-6);
    imagesc(t1,f1(f1<winLim(2) & f1>winLim(1)),s4_1);
    set(gca,'YDir','normal')
    xlabel('Time, s')
    ylabel('Frequency, Hz')
    title('Amplitude spectrogram of Sum(1-3)')
    handl = colorbar;
    colormap(jet)
    ylabel(handl, 'Magnitude, dB')
end
%% Feature Extraction for Classification:
close all;
cont = [];
showGraphs = true;
signalDetected = false;
wPlus = 250;        %-% Value by which to increase window length
winJump = 125;      %-% Data points to skip after each iteration. 
maxWinL = 250;     %-% 5s max
mW = 1:winJump:(ln - maxWinL); 
% mW = 4001:winJump:(ln-maxWinL); %TEMPORARY!!!!
ftr=1;
for i=1:length(mW)
    cWSize = 250;           %-% Start with a window size of 1s
    start = mW(i);          %-% Where to start window
    fin   = (mW(i)+(cWSize));   %-% Signal ends at start+current Win Length
    fprintf('Current index = [%d to %d]\r\n',start, fin);
    chw{1} = ch1(start:fin);  %-% temporary window variable
    chw{2} = ch2(start:fin);
    chw{3} = ch3(start:fin);
    for c = 1:3
        fch(c,:) = eegcfilt(chw{c});
    end
    F(i,:) = featureExtractionSSVEPtemp(fch(1,:), fch(2,:), fch(3,:), Fs, true);
    
    if isempty(cont)
        commandwindow;
        cont = input('Approve/continue?\n');
        clf(1);
    end
    % Feature selection: 1=accept, 9=run entirely, anything else =
    % continue/reject
    if (cont==1)
        tXSSVEP(i,:) = F(i,:);
        tY(i,1) = str2double(CLASS);
        ftr = ftr + 1;
        cont = [];
    else
        tXSSVEP(i,:) = F(i,:);
        tY(i,1) = 0; %REJECT CLASS
    end
    if ~isempty(cont) && cont~=9
        cont = [];
    end
end

filename = ['tXSSVEP_' CLASS VERSION];
prompt = ['SAVE FILE: ' filename '?\n'];
commandwindow;
cont = input(prompt);
if cont == 1
    save(filename);
end
%% Training Data Control 
% use this section to load all the tXSSVEP Data and Combine into a single
% tXSSVEP and tYSSVEP. 

% tXSSVEPALL = [;];
% tYSSVEPALL = [;];
%% CCA Test: (using CCA2)
%{
close all;clc;clear all;
load carbig;
X = [Displacement Horsepower Weight Acceleration MPG];
nans = sum(isnan(X),2) > 0;
V1 = X(~nans,1:3);
V2 = X(~nans,4:5);
[A,B,r,U,V] = CCA(V1,V2);
hold on;
plot(U(:,1),V(:,1),'.')
plot(U(:,2),V(:,2),'.')
xlabel('0.0025*Disp+0.020*HP-0.000025*Wgt')
ylabel('-0.17*Accel-0.092*MPG')
%}
