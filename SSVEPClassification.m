% %% SSVEP Classification
% The puspose of this program is to differentiate between different SSVEP
% signals including signal changes. This is accomplished using a variable
% windows. 
    %Clear 
clear;close all;clc;
    %Import Data:
ChannelNames = {['Fp1' 'Fp2' 'Fpz' 'REye']};
load('mssvep_16.6_3.mat');
% load('mssvep_15_1.mat');
% load('mssvep_10_2.mat');
% load('mssvep_12.5_1.mat')
% load('mssvep_t1_baseline');
% load('mssvep_t2_16_1.mat');

%--- LOAD CLASS ---%
CLASS = 16;
%---
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

showSpect = 1;  
%%     TODO: PLOT SPECTROGRAMS USING TWO METHODS, AND COMPARE.
f2 = figure(2);
set(f2, 'Position', [100, 100, 1200, 675]);
% Spectrograms:
if showSpect == 1
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
end   

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

%% Feature Extraction for Classification:
close all;
cont = [];
showGraphs = true;
signalDetected = false;
wPlus = 250;        %-% Value by which to increase window length
winJump = 250;      %-% Data points to skip after each iteration. 
maxWinL = 1000;     %-% 5s max
mW = 1:winJump:(ln - maxWinL); 
ftr=1;
for i=1:length(mW)
    cWSize = 250;           %-% Start with a window size of 1s
    start = mW(i);          %-% Where to start window
    fin   = mW(i)+cWSize;   %-% Signal ends at start+current Win Length
    if mod(fin-start,2)==1
        fin = fin+1;
    end
    chw{1} = ch1(start:fin);  %-% temporary window variable
    chw{2} = ch2(start:fin);
    chw{3} = ch3(start:fin);
    for c = 1:3
        fch(c,:) = eegcfilt(chw{c});
    end
    F1(i,:) = featureExtractionSSVEP(fch(1,:), fch(2,:), fch(3,:), Fs, true);
    
    if isempty(cont)
        cont = input('Approve/continue?\n');
        clf(1);
    end
    % Feature selection: 
    if (cont==1)
        F(ftr,:) = F1(i,:);
        ftr = ftr + 1;
        cont = [];
    end
end

%% Analysis *(finished, converted to function).
close all;
cont = [];
showGraphs = true;
signalDetected = false;
wPlus = 250;        %-% Value by which to increase window length
winJump = 250;      %-% Data points to skip after each iteration. 
maxWinL = 1000;     %-% 5s max
wLFFT(1,:) = [9.6 10.4];%-% windows around certain target frequencies
wLFFT(2,:) = [11.9 12.7]; 
wLFFT(3,:) = [14.6 15.5];
wLFFT(4,:) = [16.2 16.7];
wLPSD(1,:) = [9.9 10.1];
wLPSD(2,:) = [12 13]; 
wLPSD(3,:) = [14.9 15.1];
wLPSD(4,:) = [16 17];
winLim = [9 17.6];
xL = [9.0 17.2];    %-% Xlimit to use for plotting
yL = [8 18];
fL = [8.0 17.9];    %-% Frequency range to filter
oN = 7;             %-% Filter Order (BW)
mW = 1:winJump:(ln - maxWinL);   %-% Separate into moving Windows (mW)
cWSize = 250;       %-% Current Window Size
if isempty(cont)    %-% Set up plot
    fH = figure(1); %-% Figure Handle
    set(fH, 'Position', [0, 0, 1280, 920]);
%     fH2 = figure(2);
%     set(fH2, 'Position', [960, 0, 960, 920]);
end
nCh = 3;            %-% For now will use Fp1, Fp2 and Fpz
chw = cell(nCh,1);
fchw = cell(nCh,1);
Ch = cell(nCh,1);
% f = cell(length(mW), 1);
clear f;
clear f2;
fPSD = cell(length(mW), 1);
str{1} = '^c';
str{2} = '^r';
str{3} = '^m';
 %%% Create different boolean conditions that have to be enabled. 
ftr = 1;
for i=1:length(mW)
    cWSize = 250;           %-% Start with a window size of 1s
    start = mW(i);          %-% Where to start window
    fin   = mW(i)+cWSize;   %-% Signal ends at start+current Win Length
    if mod(fin-start,2)==1
        fin = fin+1;
    end
    chw{1} = ch1(start:fin);  %-% temporary window variable
    chw{2} = ch2(start:fin);
    chw{3} = ch3(start:fin);
    %Filter Normally:
    for ch=1:nCh        %-% Channels 1-3 or however many we use.
        fchw{ch} = customFilt(chw{ch}, Fs, fL, oN);
        % #1 Take FFT:
        [f, Ch{ch}.FFT{i}] = get_nfft_data(fchw{ch}, Fs, 2048);
            % #1.1 Find Peaks and M/I
            [Ch{ch}.FFT_PKS{i}, Ch{ch}.FFT_L{i}] = findpeaks(Ch{ch}.FFT{i},'SortStr','descend');
            Ch{ch}.FFT_L{i} = Ch{ch}.FFT_L{i}(:);
            if length(Ch{ch}.FFT_PKS{i})>1
                %Peak max minus min
                Ch{ch}.FFT_Ltop{i} = f(Ch{ch}.FFT_L{i}(1:2,1));
                for w = 1:4
                    if Ch{ch}.FFT_Ltop{i}(1,1)>wLFFT(w,1) && Ch{ch}.FFT_Ltop{i}(1,1)<wLFFT(w,2) 
                        Ch{ch}.FFT_MMM{i} = Ch{ch}.FFT_PKS{i}(1) - Ch{ch}.FFT_PKS{i}(2);
                        Ch{ch}.FFT_PkRatio{i} = Ch{ch}.FFT_PKS{i}(1)/Ch{ch}.FFT_PKS{i}(2);
                        Ch{ch}.wLFFT{i} = w;
                        break;
                    else
                        Ch{ch}.FFT_MMM{i} = 0;
                        Ch{ch}.FFT_PkRatio{i} = 0;
                        Ch{ch}.wLFFT{i} = 0;
                    end
                end
            end
        % #2 Take PSD Estimate: (Welch method)
        [Ch{ch}.PSD{i}, fPSD] = welch_psd(fchw{ch}, Fs, hann(fin-start));%fin-start
            % #2.2 Find Peaks and Max
            [Ch{ch}.PSD_PKS{i}, Ch{ch}.PSD_L{i}] = findpeaks(Ch{ch}.PSD{i},'SortStr','descend');
            Ch{ch}.PSD_L{i} = Ch{ch}.PSD_L{i}(:);
            if length(Ch{ch}.PSD_PKS{i})>1
                Ch{ch}.PSD_Ltop{i} = fPSD(Ch{ch}.PSD_L{i}(1:2,1));
                for w = 1:4
                    if Ch{ch}.PSD_Ltop{i}(1,1)>=wLPSD(w,1) && Ch{ch}.PSD_Ltop{i}(1,1)<=wLPSD(w,2)
                        Ch{ch}.PSD_MMM{i} = Ch{ch}.PSD_PKS{i}(1) - Ch{ch}.PSD_PKS{i}(2);
                        Ch{ch}.PSD_PkRatio{i} = Ch{ch}.PSD_PKS{i}(1) / Ch{ch}.PSD_PKS{i}(2);
                        Ch{ch}.wLPSD{i} = w;
                        break;
                    else
                        Ch{ch}.PSD_MMM{i} = 0;
                        Ch{ch}.PSD_PkRatio{i} = 0;
                        Ch{ch}.wLPSD{i} = 0;
                    end
                end
            end
        if(showGraphs)
            figure(1)
            subplot(2,2,1)
            hold on;
            plot(f,Ch{ch}.FFT{i}),xlim(xL);
            %Pks and max:
            plot(f(Ch{ch}.FFT_L{i}), Ch{ch}.FFT_PKS{i},str{ch});
            title('FFT');
            subplot(2,2,2)
            hold on;
            plot(fPSD, Ch{ch}.PSD{i}),xlim(xL);
            plot(fPSD(Ch{ch}.PSD_L{i}), Ch{ch}.PSD_PKS{i},str{ch});
            title('PSD');
        end
    end
    Ch{4}.FFT{i} = (Ch{1}.FFT{i}+Ch{2}.FFT{i}+Ch{3}.FFT{i});
    [Ch{4}.FFT_PKS{i}, Ch{4}.FFT_L{i}] = findpeaks(Ch{4}.FFT{i},'SortStr','descend');
    if length(Ch{4}.FFT_PKS{i})>1
        Ch{4}.FFT_L{i} = Ch{4}.FFT_L{i}(:);
        Ch{4}.FFT_Ltop{i} = f(Ch{4}.FFT_L{i}(1:2,1));
        for w = 1:4
            if Ch{4}.FFT_Ltop{i}(1,1)>wLFFT(w,1) && Ch{4}.FFT_Ltop{i}(1,1)<wLFFT(w,2) 
                Ch{4}.FFT_MMM{i} = Ch{4}.FFT_PKS{i}(1) - Ch{4}.FFT_PKS{i}(2);
                Ch{4}.FFT_PkRatio{i} = Ch{4}.FFT_PKS{i}(1)/Ch{4}.FFT_PKS{i}(2);
                Ch{4}.wLFFT{i} = w;
                break;
            else
                Ch{4}.FFT_MMM{i} = 0;
                Ch{4}.FFT_PkRatio{i} = 0;
                Ch{4}.wLFFT{i} = 0;
            end
        end
    end
    Ch{4}.PSD{i} = Ch{1}.PSD{i}+Ch{2}.PSD{i}+Ch{3}.PSD{i};
    [Ch{4}.PSD_PKS{i}, Ch{4}.PSD_L{i}] = findpeaks(Ch{4}.PSD{i},'SortStr','descend');
    if length(Ch{4}.PSD_PKS{i})>1
        Ch{4}.PSD_L{i} = Ch{4}.PSD_L{i}(:);
        Ch{4}.PSD_Ltop{i} = fPSD(Ch{4}.PSD_L{i}(1:2,1));
        for w = 1:4
            if Ch{4}.PSD_Ltop{i}(1,1)>=wLPSD(w,1) && Ch{4}.PSD_Ltop{i}(1,1)<=wLPSD(w,2)
                Ch{4}.PSD_MMM{i} = Ch{4}.PSD_PKS{i}(1) - Ch{4}.PSD_PKS{i}(2);
                Ch{4}.PSD_PkRatio{i} = Ch{4}.PSD_PKS{i}(1) / Ch{4}.PSD_PKS{i}(2);
                Ch{4}.wLPSD{i} = w;
                break;
            else
                Ch{4}.PSD_MMM{i} = 0;
                Ch{4}.PSD_PkRatio{i} = 0;
                Ch{4}.wLPSD{i} = 0;
            end
        end
    end
    %Alignment:
    %averageFFTPeak is for KNN classification. 
    % b1 = all maxfreqs (FFT) for this iteration triggered a known frequency. 
    for chn = 1:4
        FFTPeaks1(chn, i) = Ch{chn}.FFT_Ltop{1,i}(1);
        FFTPeaks2(chn, i) = Ch{chn}.FFT_Ltop{1,i}(2);
    end
    averageFFTPeak{i} = mean([Ch{1}.FFT_Ltop{1,i}(1) Ch{2}.FFT_Ltop{1,i}(1) ...
        Ch{3}.FFT_Ltop{1,i}(1) Ch{4}.FFT_Ltop{1,i}(1)]);
    averageFFTPeak2{i} = mean([Ch{1}.FFT_Ltop{1,i}(2) Ch{2}.FFT_Ltop{1,i}(2) ...
        Ch{3}.FFT_Ltop{1,i}(2) Ch{4}.FFT_Ltop{1,i}(2)]);
    b1(i) = (Ch{1}.wLFFT{1,i}~=0) && (Ch{2}.wLFFT{1,i}~=0) && (Ch{3}.wLFFT{1,i}~=0) && (Ch{4}.wLFFT{1,i}~=0);
    averagePSDPeak{i} = mean([Ch{1}.PSD_Ltop{1,i}(1) Ch{2}.PSD_Ltop{1,i}(1) ...
        Ch{3}.PSD_Ltop{1,i}(1) Ch{4}.PSD_Ltop{1,i}(1)]);
    b2(i) = (Ch{1}.wLPSD{1,i}~=0) && (Ch{2}.wLPSD{1,i}~=0) && (Ch{3}.wLPSD{1,i}~=0) && (Ch{4}.wLPSD{1,i}~=0);
%     if ~b1{i}
%         check for tolerance (either 1 or two
%     end
    if(showGraphs)
        figure(1)
        hold off;
        subplot(2,2,3);hold on;        
        plot(f,Ch{4}.FFT{i}),xlim(xL);
        plot(f(Ch{4}.FFT_L{i}),Ch{4}.FFT_PKS{i},'^r');
        hold off;
        title('FFT (P1+P2+P3)');
        subplot(2,2,4);hold on;
        plot(fPSD,Ch{4}.PSD{i}),xlim(xL);
        plot(fPSD(Ch{4}.PSD_L{i}),Ch{4}.PSD_PKS{i},'^r');
        title('PSD (P1+P2+P3)');
        hold off;
        hold off;
    end
    %%% TODO: analyse individual windows:
    %%% SET CONDITIONAL STATEMENTS TO OUTPUT [10 12 15 16], if NOT
    %%% DETECTABLE, SKIP AND INCREASE WIN SIZE
    if isempty(cont)
        cont = input('Approve/continue?\n');
        clf(fH);
    end
    % Feature selection: 
    if (cont==1)
        %{N, 1(250samples)}
%         ft_ch1 = 
%         ft_ch2 = 
%         ft_ch3 = 
%         ft_ch4 = 
%         F{ftr,1} = [ft_ch1, ft_ch2, ft_ch3, ft_ch4];
        ftr = ftr + 1;
        cont = [];
    end
    %%% TODO: Save final win length
    % Cannot find match
    if ~signalDetected % -- Temporary & Redundant
        cWSize = cWSize + wPlus;
    end
end
%% Part 2: 500 & longer
%{
%}
%-% STFT Variables:
wlen = 256;         %-% Length of the hamming window
h = 64;             %-% hop size
% win = hamming(wlen, 'periodic'); %CAN BE CONVERTED TO C IF WLEN IS CONSTANT
nfft = 2048;
for i=1:length(mW)
    cWSize = 250;           %-% Start with a window size of 1s
    while cWSize<maxWinL %~signalDetected   %-% Continue iterating until signal has been detected. 
        cWSize
        start = mW(i);      %-% Where to start window
        fin   = mW(i)+cWSize;   %-% Signal ends at start+current Win Length
        if mod(fin-start,2)==1
            fin = fin+1;
        end
        chw{1} = ch1(start:fin);  %-% temporary window variable
        chw{2} = ch2(start:fin);
        chw{3} = ch3(start:fin);
        %Filter Normally:
        for ch=1:nCh        %-% Channels 1-3 or however many we use.
            fchw{ch} = customFilt(chw{ch}, Fs, fL, oN);
            % #1 Take FFT:
            [f{i}, Ch{ch}.FFT{i}] = get_fft_data(fchw{ch}, Fs);
                % #1.1 Find Peaks and M/I
                [Ch{ch}.FFT_MAX{i}, Ch{ch}.FFT_I{i}] = max(Ch{ch}.FFT{i});
                [Ch{ch}.FFT_PKS{i}, Ch{ch}.FFT_L{i}] = findpeaks(Ch{ch}.FFT{i});
            % #2 Take PSD Estimate: (Welch method)
            [Ch{ch}.PSD{i}, fPSD{i}] = welch_psd(fchw{ch}, Fs, hann(fin-start));
                % #2.2 Find Peaks and Max
                [Ch{ch}.PSD_MAX{i}, Ch{ch}.PSD_I{i}] = max(Ch{ch}.PSD{i});
                [Ch{ch}.PSD_PKS{i}, Ch{ch}.PSD_L{i}] = findpeaks(Ch{ch}.PSD{i});
            if(cWSize >= 500)
                % #3 Take STFT.
                K = sum(hamming(wlen, 'periodic'))/wlen;
                [Ch{ch}.sSTFT{i}, Ch{ch}.fSTFT{i}, Ch{ch}.tSTFT{i}] = stft(fchw{ch},h,nfft,Fs);
                Ch{ch}.sSTFTlog{i} = 20*log10(abs(Ch{ch}.sSTFT{i})/wlen/K + 1E-6);
                % Separate into windows (boxed around target frequencies).
            end
            if isempty(cont)
                figure(1)
                subplot(4,2,5)
                hold on;
                plot(f{i},Ch{ch}.FFT{i}),xlim(xL);
                %Pks and max:
                plot(f{i}(Ch{ch}.FFT_I{i}), Ch{ch}.FFT_MAX{i},'or');
                plot(f{i}(Ch{ch}.FFT_L{i}), Ch{ch}.FFT_PKS{i},str{ch});
                title('FFT');
                subplot(4,2,6)
                hold on;
                plot(fPSD{i}, Ch{ch}.PSD{i}),xlim(xL);
                plot(fPSD{i}(Ch{ch}.PSD_I{i}), Ch{ch}.PSD_MAX{i},'or');
                plot(fPSD{i}(Ch{ch}.PSD_L{i}), Ch{ch}.PSD_PKS{i},str{ch});
                title('PSD');
                if(cWSize >= 500)
                    subplot(4,2,ch)
                    imagesc(Ch{ch}.tSTFT{i},Ch{ch}.fSTFT{i},Ch{ch}.sSTFTlog{i}),ylim(yL);
                    set(gca,'YDir','normal')
                    xlabel('Time, s')
                    ylabel('Frequency, Hz')
                    title('Amplitude spectrogram of the signal')
                    handl = colorbar;
                    colormap(jet);
                    ylabel(handl, 'Magnitude, dB')
                end
            end
        end
        if (cWSize>=500)
            ChCombined.sSTFT{i,1} = 20*log10(abs(Ch{1}.sSTFT{i}(Ch{1}.fSTFT{i}<winLim(2) & Ch{1}.fSTFT{i}>winLim(1),:))/wlen/K + 1E-6)+...
                        20*log10(abs(Ch{2}.sSTFT{i}(Ch{2}.fSTFT{i}<winLim(2) & Ch{2}.fSTFT{i}>winLim(1),:))/wlen/K + 1E-6)+...
                        20*log10(abs(Ch{3}.sSTFT{i}(Ch{3}.fSTFT{i}<winLim(2) & Ch{3}.fSTFT{i}>winLim(1),:))/wlen/K + 1E-6);
            ChCombined.SummedRows{i} = scaleAbs(sum(ChCombined.sSTFT{i},2));

            figure(2)
            subplot(2,1,1);
            hold on;
            ChCombined.fSTFT{i} = Ch{1}.fSTFT{i}(Ch{1}.fSTFT{i}<winLim(2) & Ch{1}.fSTFT{i}>winLim(1));
            plot(ChCombined.fSTFT{i}, ChCombined.SummedRows{i});
            [ChCombined.Max{i}, ChCombined.MaxInd{i}] = max( ChCombined.SummedRows{i} );
            plot(ChCombined.fSTFT{i}(ChCombined.MaxInd{i}),ChCombined.Max{i},'-or');
            hold off;
        end
        if isempty(cont)
            figure(1)
            hold off;
            subplot(4,2,7);
            Ch{4}.FFT{i} = (Ch{1}.FFT{i}+Ch{2}.FFT{i}+Ch{3}.FFT{i});
            plot(f{i},Ch{4}.FFT{i}),xlim(xL);
            title('FFT (P1+P2+P3)');
            subplot(4,2,8);
            plot(fPSD{i},(Ch{1}.PSD{i}+Ch{2}.PSD{i}+Ch{3}.PSD{i})),xlim(xL);
            title('PSD (P1+P2+P3)');
            hold off;
            if(cWSize >= 500)
                subplot(4,2,4);
                imagesc(Ch{1}.tSTFT{i},Ch{1}.fSTFT{i}(Ch{1}.fSTFT{i}<winLim(2) & Ch{1}.fSTFT{i}>winLim(1)),ChCombined.sSTFT{i}),ylim(winLim);
                set(gca,'YDir','normal')
                xlabel('Time, s')
                ylabel('Frequency, Hz')
                title('Amplitude spectrogram of combined signal')
                handl = colorbar;
                ylabel(handl, 'Magnitude, dB')
            end
        end
        %%% TODO: analyse individual windows:
        %%% SET CONDITIONAL STATEMENTS TO OUTPUT [10 12 15 16], if NOT
        %%% DETECTABLE, SKIP AND INCREASE WIN SIZE
        if isempty(cont)
            cont = input('Approve/continue?\n');
            clf(fH);
        end
        %%% TODO: Save final win length
        % Cannot find match
        if ~signalDetected % -- Temporary & Redundant
            cWSize = cWSize + wPlus;
        end
    end
end

%% Analysis (OLD)
% start with smallest possible window:
% TODO: EACH winL will result in a feature. USE: FFT, STFT, PSD, and ??
% SEPARATE FEATURES BY 1s, 2s, and 4s elapsed and have separate classifiers
% for each. 
cccc = input('Run Next Section?? (Ctrl-C to cancel)\n');
cont = [];
close all;
minlen = min([ length(ch1) length(ch2) length(ch3) ]);
winL = [ 250 376 500 626 750 876 1000 ]; %0.5?4s
newWin = 250;
ii = 1:newWin:(minlen-max(winL));
xl = [5 25];
wlen = 2^nextpow2(Fs);
h=wlen/4;
win = hamming(wlen, 'periodic'); %CAN BE CONVERTED TO C IF WLEN IS CONSTANT
nfft = 2^nextpow2(wlen+1);
K = sum(hamming(wlen, 'periodic'))/wlen;
if isempty(cont)
    fH = figure(1);
    set(fH, 'Position', [100, 100, 1600, 1000]);
end
loc500 = find(winL==500)-1;
% ----- CLASSIFIER THRESHOLDS ----- %
PeakThresholdPSD = 0.5E-11;
% -- Actual Frequencies: (rounded to 4 decimal places)
f0 = [10.0000 12.5000 15.1515 16.6667];
% Tolerance (f0 +/- f0t)
% Low-res tolerances.
f0t_low = 0.45;
f0r_low = [f0'-f0t_low, f0'+f0t_low];
% High-res tolerances.
f0t = 0.3;
f0r = [f0'-f0t, f0'+f0t];
% Frequency ranges (based on tolerance, for STFT.):
% -- Classifier Possible Outputs:
fp = [10 12 15 16];
tYrejectCount = 1;
for i = 1:length(ii)
    for j = 1:length(winL)
        if isempty(cont)
            fprintf('%d -> %d \n', ii(i),ii(i)+winL(j)-1);
        end
        %TODO: REPLACE CUSTOM FILTER WITH STATIC ONE FOR CONVERSION
        Ch1.Windows{i,j} = customFilt( ch1(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [Ch1.fFFT{j}, Ch1.FFT{i,j}] = get_fft_data(Ch1.Windows{i,j}, Fs);
        [Ch1.MaxFFT{i,j}, Ch1.IndicesMaxFFT{i,j}] = max(Ch1.FFT{i,j});
        
        Ch2.Windows{i,j} = customFilt( ch2(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [~, Ch2.FFT{i,j}] = get_fft_data(Ch2.Windows{i,j}, Fs);
        [Ch2.MaxFFT{i,j}, Ch2.IndicesMaxFFT{i,j}] = max(Ch2.FFT{i,j});
        
        Ch3.Windows{i,j} = customFilt( ch3(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [~, Ch3.FFT{i,j}] = get_fft_data(Ch3.Windows{i,j}, Fs);
        [Ch3.MaxFFT{i,j}, Ch3.IndicesMaxFFT{i,j}] = max(Ch3.FFT{i,j});
            %%% --- POWER SPECTRAL DENSITY EST --- %%%
                % Note: DOES NOT ACCEPT WINDOWS OF ODD LENGTH %
            [Ch1.PSDData{i,j}, Ch1.fPSD{j}] = welch_psd(Ch1.Windows{i,j}, Fs, hann(winL(j)));
            [Ch2.PSDData{i,j}, Ch2.fPSD{j}] = welch_psd(Ch2.Windows{i,j}, Fs, hann(winL(j)));
            [Ch3.PSDData{i,j}, Ch3.fPSD{j}] = welch_psd(Ch3.Windows{i,j}, Fs, hann(winL(j)));
                % FIND MAX VALUE:
            [Ch1.PSDPeak{i,j}, Ch1.PSDi{i,j}] = max(Ch1.PSDData{i,j});
            [Ch2.PSDPeak{i,j}, Ch2.PSDi{i,j}] = max(Ch2.PSDData{i,j});
            [Ch3.PSDPeak{i,j}, Ch3.PSDi{i,j}] = max(Ch3.PSDData{i,j});
                % FIND LOCAL MAX VALUES:
                    % IGNORE EVERYTHING OUTSIDE OF [9 18]Hz
                    % TODO: FIND LOCAL MAX IN FOUR FREQ REGIONS:
            if isempty(cont)
                %%% -------- PLOT PSDs ------------ %%%
                subplot(3,3,[4 6]); hold on;
                plot(Ch1.fPSD{j}, Ch1.PSDData{i,j}),xlim(xl);
                plot(Ch2.fPSD{j}, Ch2.PSDData{i,j}),xlim(xl);
                plot(Ch3.fPSD{j}, Ch3.PSDData{i,j}),xlim(xl);
                % Plot Max Values
                plot(Ch1.fPSD{j}(Ch1.PSDi{i,j}),Ch1.PSDPeak{i,j},'or');
                plot(Ch2.fPSD{j}(Ch2.PSDi{i,j}),Ch2.PSDPeak{i,j},'om');
                plot(Ch3.fPSD{j}(Ch3.PSDi{i,j}),Ch3.PSDPeak{i,j},'og');
                str = [' f = ' num2str(Ch1.fPSD{j}(Ch1.PSDi{i,j}))];
                text(Ch1.fPSD{j}(Ch1.PSDi{i,j}),Ch1.PSDPeak{i,j},str);
                str = [' f = ' num2str(Ch2.fPSD{j}(Ch2.PSDi{i,j}))];
                text(Ch2.fPSD{j}(Ch2.PSDi{i,j}),Ch2.PSDPeak{i,j},str);
                str = [' f = ' num2str(Ch3.fPSD{j}(Ch3.PSDi{i,j}))];
                text(Ch3.fPSD{j}(Ch3.PSDi{i,j}),Ch3.PSDPeak{i,j},str);
                % Plot Local Maxima
                xlabel('Normalized frequency'); %ylabel('PSD [dB]');
                ylabel('Power Spectrum Magnitude');
                title('Power Spectral Test');
                hold off;
                %%% -------- PLOT FFTs ------------ %%%
                subplot(3,3,[1 3]);
                hold on;
                plot(Ch1.fFFT{j}, Ch1.FFT{i,j}),xlim(xl);
                    plot(Ch1.fFFT{j}(Ch1.IndicesMaxFFT{i,j}), Ch1.MaxFFT{i,j},'-.r*');
                    str = [' f = ' num2str(Ch1.fFFT{j}(Ch1.IndicesMaxFFT{i,j})) '  M = ' num2str( Ch1.MaxFFT{i,j} )];
                    text(Ch1.fFFT{j}(Ch1.IndicesMaxFFT{i,j}), Ch1.MaxFFT{i,j}, str);
                plot(Ch1.fFFT{j}, Ch2.FFT{i,j}),xlim(xl);
                    plot(Ch1.fFFT{j}(Ch2.IndicesMaxFFT{i,j}), Ch2.MaxFFT{i,j},'-.m*');
                    str = [' f = ' num2str(Ch1.fFFT{j}(Ch2.IndicesMaxFFT{i,j})) '  M = ' num2str( Ch2.MaxFFT{i,j} )];
                    text(Ch1.fFFT{j}(Ch2.IndicesMaxFFT{i,j}), Ch2.MaxFFT{i,j}, str);
                plot(Ch1.fFFT{j}, Ch3.FFT{i,j}),xlim(xl);
                    plot(Ch1.fFFT{j}(Ch3.IndicesMaxFFT{i,j}), Ch3.MaxFFT{i,j},'-.c*');
                    str = [' f = ' num2str(Ch1.fFFT{j}(Ch3.IndicesMaxFFT{i,j})) '  M = ' num2str( Ch3.MaxFFT{i,j} )];
                    text(Ch1.fFFT{j}(Ch3.IndicesMaxFFT{i,j}), Ch3.MaxFFT{i,j}, str);
                title('FFT (Ch 1-3): With Peaks');
                ylabel('|P1(f)|');
                xlabel('f (Hz)');
                hold off;
            end
            %%% --- APPLY STFT --- %%%
            if(winL(j) >= 500)
                %CH1
                [Ch1.sSTFT{i, j-loc500}, Ch1.fSTFT{j-loc500}, Ch1.tSTFT{j-loc500}] = stft( Ch1.Windows{i,j}, h, nfft, Fs );
                Ch1.sSTFT{i, j-loc500} = 20*log10(abs(Ch1.sSTFT{i, j-loc500})/wlen/K + 1e-6); 
                %CH2
                [Ch2.sSTFT{i, j-loc500}, Ch2.fSTFT{j-loc500}, Ch2.tSTFT{j-loc500}] = stft( Ch2.Windows{i,j}, h, nfft, Fs );
                Ch2.sSTFT{i, j-loc500} = 20*log10(abs(Ch2.sSTFT{i, j-loc500})/wlen/K + 1e-6);
                %CH3
                [Ch3.sSTFT{i, j-loc500}, Ch3.fSTFT{j-loc500}, Ch3.tSTFT{j-loc500}] = stft( Ch3.Windows{i,j}, h, nfft, Fs );
                Ch3.sSTFT{i, j-loc500} = 20*log10(abs(Ch3.sSTFT{i, j-loc500})/wlen/K + 1e-6);
                %%%%%%-TODO feature extraction from s f t [5->20Hz]
                if isempty(cont)
                    subplot(3,3,7);
                    imagesc(Ch1.tSTFT{j-loc500},Ch1.fSTFT{j-loc500},Ch1.sSTFT{i, j-loc500}),ylim([7.25 20]);
                    set(gca,'YDir','normal')
                    xlabel('Time, s')
                    ylabel('Frequency, Hz')
                    title('Amplitude spectrogram of the signal')
                    handl = colorbar;
                    ylabel(handl, 'Magnitude, dB')
                    subplot(3,3,8);
                    imagesc(Ch2.tSTFT{j-loc500},Ch2.fSTFT{j-loc500},Ch2.sSTFT{i, j-loc500}),ylim([7.25 20]);
                    set(gca,'YDir','normal')
                    xlabel('Time, s')
                    ylabel('Frequency, Hz')
                    title('Amplitude spectrogram of the signal')
                    handl = colorbar;
                    ylabel(handl, 'Magnitude, dB')
                    subplot(3,3,9);
                    imagesc(Ch3.tSTFT{j-loc500},Ch3.fSTFT{j-loc500},Ch3.sSTFT{i, j-loc500}),ylim([7.25 20]);
                    set(gca,'YDir','normal')
                    xlabel('Time, s')
                    ylabel('Frequency, Hz')
                    title('Amplitude spectrogram of the signal')
                    handl = colorbar;
                    ylabel(handl, 'Magnitude, dB')
                end
            end
        % --- DECISION TREE --- %:
        %1 - Boolean Checks:
        for b = 1:length(f0)
            if winL(j)<=376
                frange = f0r_low(b,:);
                Ch1.B{i,j}.b1(b) = isWithin(Ch1.fPSD{j}(Ch1.PSDi{i,j}),frange);
                Ch2.B{i,j}.b1(b) = isWithin(Ch2.fPSD{j}(Ch1.PSDi{i,j}),frange);
                Ch3.B{i,j}.b1(b) = isWithin(Ch3.fPSD{j}(Ch1.PSDi{i,j}),frange);
                if(Ch1.B{i,j}.b1(b))
                    Ch1.B{i,j}.m1(b) = Ch1.PSDPeak{i,j};
                else
                    Ch1.B{i,j}.m1(b) = 0.0;
                end
                if(Ch2.B{i,j}.b1(b))
                    Ch2.B{i,j}.m1(b) = Ch2.PSDPeak{i,j};
                else
                    Ch2.B{i,j}.m1(b) = 0.0;
                end
                if(Ch3.B{i,j}.b1(b))
                    Ch3.B{i,j}.m1(b) = Ch3.PSDPeak{i,j};
                else
                    Ch3.B{i,j}.m1(b) = 0.0;
                end
            else
                frange = f0r(b,:);
                Ch1.B{i,j}.b1(b) = isWithin(Ch1.fPSD{j}(Ch1.PSDi{i,j}),frange);
                Ch2.B{i,j}.b1(b) = isWithin(Ch2.fPSD{j}(Ch1.PSDi{i,j}),frange);
                Ch3.B{i,j}.b1(b) = isWithin(Ch3.fPSD{j}(Ch1.PSDi{i,j}),frange);
                if(Ch1.B{i,j}.b1(b))
                    Ch1.B{i,j}.m1(b) = Ch1.PSDPeak{i,j};
                else
                    Ch1.B{i,j}.m1(b) = 0.0;
                end
                if(Ch2.B{i,j}.b1(b))
                    Ch2.B{i,j}.m1(b) = Ch2.PSDPeak{i,j};
                else
                    Ch2.B{i,j}.m1(b) = 0.0;
                end
                if(Ch3.B{i,j}.b1(b))
                    Ch3.B{i,j}.m1(b) = Ch3.PSDPeak{i,j};
                else
                    Ch3.B{i,j}.m1(b) = 0.0;
                end
            end
        end
        % --- USER INPUT --- %
        if isempty(cont)
            cont = input('Approve/continue?\n');
            if ~isempty(cont)
                if cont~=0
                    tYreject(tYrejectCount,:) = i;
                    tYrejectCount = tYrejectCount+1;
                    cont = [];
                end
            end
            clf(1);
        end
    end
end

%% Combine Features:
%Preallocate:
tY = cell(length(winL),1);
tXCh1 = cell(length(winL),1);
tXCh2 = tXCh1;
tXCh3 = tXCh1;
tX = tXCh1;
for w = 1:length(winL)
    for r = 1:size(Ch1.MaxFFT,1) %Rows {w}
%         Ch1:
        tXCh1{w}.FFT(r,1) = Ch1.fFFT{w}(Ch1.IndicesMaxFFT{r,w});
        tXCh1{w}.FFT(r,2) = Ch1.MaxFFT{r,w};
        tXCh1{w}.PSD(r,1) = Ch1.fPSD{w}(Ch1.PSDi{r,w});
        tXCh1{w}.PSD(r,2) = Ch1.PSDPeak{r,w};
%         Ch2:
        tXCh2{w}.FFT(r,1) = Ch1.fFFT{w}(Ch2.IndicesMaxFFT{r,w});
        tXCh2{w}.FFT(r,2) = Ch2.MaxFFT{r,w};
        tXCh2{w}.PSD(r,1) = Ch2.fPSD{w}(Ch2.PSDi{r,w});
        tXCh2{w}.PSD(r,2) = Ch2.PSDPeak{r,w};
%         Ch3:
        tXCh3{w}.FFT(r,1) = Ch1.fFFT{w}(Ch3.IndicesMaxFFT{r,w});
        tXCh3{w}.FFT(r,2) = Ch3.MaxFFT{r,w};
        tXCh3{w}.PSD(r,1) = Ch3.fPSD{w}(Ch3.PSDi{r,w});
        tXCh3{w}.PSD(r,2) = Ch3.PSDPeak{r,w};
        for c = 1:length(f0)
            tXCh1{w}.M(r, c) = Ch1.B{r,w}.m1(c);
            tXCh2{w}.M(r, c) = Ch2.B{r,w}.m1(c);
            tXCh3{w}.M(r, c) = Ch3.B{r,w}.m1(c);
        end
        tY{w}(r,1) = CLASS; 
    end
    tX{w} = [tXCh1{w}.FFT tXCh1{w}.PSD tXCh1{w}.M ...
    tXCh2{w}.FFT tXCh2{w}.PSD tXCh2{w}.M ...
    tXCh3{w}.FFT tXCh3{w}.PSD tXCh3{w}.M ];
end
tXtY = [tX tY];
numberSamplesX = length(ii);
clearvars -except tX tY tXtY winL f0 f0r f0r_low tXCh1 tXCh2 tXCh3 CLASS numberSamplesX seconds

%% Compile all Data:
clear all;clc;close all;
    % Baseline
load('tD0_baseline1.mat')
tXtY0 = tXtY;
load('tD0_baseline2.mat')
tXtY1 = tXtY;
    % 10Hz
load('tD10_1.mat')
tXtY10_1 = tXtY;
load('tD10_2.mat')
tXtY10_2 = tXtY;
    %12Hz
load('tD12_1.mat')
tXtY12_1 = tXtY;
load('tD12_2.mat')
tXtY12_2 = tXtY;
load('tD12_3.mat')
tXtY12_3 = tXtY;
    %15Hz
load('tD15_1.mat')
tXtY15_1 = tXtY;
load('tD15_2.mat')
tXtY15_2 = tXtY;
load('tD15_3.mat')
tXtY15_3 = tXtY;
    %16Hz
load('tD16.mat')
tXtY16 = tXtY;
for i = 1:size(tXtY,1)
    tXtY_C{i,1} = [tXtY0{i,1};tXtY1{i,1};...
        tXtY10_1{i,1};tXtY10_2{i,1};...
        tXtY12_1{i,1};tXtY12_2{i,1};tXtY12_3{i,1};...
        tXtY15_1{i,1};tXtY15_2{i,1};tXtY15_3{i,1};tXtY16{i,1}]; 
    tXtY_C{i,2} = [tXtY0{i,2};tXtY1{i,2};...
        tXtY10_1{i,2};tXtY10_2{i,2};...
        tXtY12_1{i,2};tXtY12_2{i,2};tXtY12_3{i,2};...
        tXtY15_1{i,2};tXtY15_2{i,2};tXtY15_3{i,2};tXtY16{i,2}];
end
clearvars -except tXtY*
clear tXtY
tXtY_250 = [tXtY_C{1,1} tXtY_C{1,2}];
tXtY_376 = [tXtY_C{2,1} tXtY_C{2,2}];
tXtY_500 = [tXtY_C{3,1} tXtY_C{3,2}];
tXtY_626 = [tXtY_C{4,1} tXtY_C{4,2}];
tXtY_750 = [tXtY_C{5,1} tXtY_C{5,2}];
tXtY_876 = [tXtY_C{6,1} tXtY_C{6,2}];
tXtY_1000 = [tXtY_C{7,1} tXtY_C{7,2}];
[Wx, Wy] = cca(tXtY_C{7,1}', tXtY_C{7,2}');
%% Sort out what will be passed to CCA. 
% etx = floor(size(tX,1)/2);
% [Wx, Wy, r] = cca(tX(1:etx,:), tX(etx+1:end,:));



 
%- END
