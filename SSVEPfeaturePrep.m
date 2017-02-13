%% Clear & Load Data
clear;close all;clc;
which_pc = input('WHICH PC? : 1=home, 0=work \n');
% load the BioRadio API using a MATLAB's .NET interface
if which_pc == 1
    dataRootFolder = 'C:\Users\Musa Mahmood\Dropbox\Public\_VCU\Yeo Lab\_SSVEP\_MATLAB-SSVEP-Classification\data';
else
    dataRootFolder = 'C:\Users\mahmoodms\Dropbox\Public\_VCU\Yeo Lab\_SSVEP\_MATLAB-SSVEP-Classification\data';
end
folder{1} = '\EOG_snap\'; 
folder{2} = '\SSVEP_snap\'; 
folder{3} = '\SSVEP_alt_setups\';
% filename = 'Dad_X1_10Hz.mat';
%     disp(filename);
% load([dataRootFolder folder{3} filename]);
load('Marc_nonHair_17.24Hz_4.mat');
% expectedFreq = input('what frequency?\n');
fp1 = Trial{1}(1:end-250,1); %ignore last second
fp2 = Trial{2}(1:end-250,1);
Fs = SamplingRate;
%%% Plot FFT (Raw)
f = [7 30];
N = 5;
fp1_f = eeg_h_custom(fp1, Fs, f, N);
fp2_f = eeg_h_custom(fp2, Fs, f, N);
[f, P1] = get_fft_data(fp1_f, Fs);
[f2, P2] = get_fft_data(fp2_f, Fs);
figure; hold on;
plot(f,P1),xlim([1 35]);
plot(f2,P2),xlim([1 35]);
hold off;
Lseconds = length(fp1_f);
h=1/250;
t=0:h:length(fp1)/250-h;
markers = trainingData{1};
    % Spectrogram:
    figure;
[~, Fspect, T, P] = spectrogram(fp1_f, 5*Fs,4*Fs,10*Fs,Fs);
imagesc(T, Fspect(Fspect<50 & Fspect>1), 10*log10(P(Fspect<50 & Fspect>1,:)));
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
cb = colorbar;
ylabel(cb, 'Power (db)')
colormap(jet)
title('Channel 1', 'FontSize', 14)
    figure;
[~, Fspect, T, P] = spectrogram(fp2_f, 5*Fs,4*Fs,10*Fs,Fs);
imagesc(T, Fspect(Fspect<50 & Fspect>1), 10*log10(P(Fspect<50 & Fspect>1,:)));
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
cb = colorbar;
ylabel(cb, 'Power (db)')
colormap(jet)
title('Channel 2', 'FontSize', 14)
%% Analysis
    %keep data selection:
close all;
%     classFreqs = [6 10 12 14];
%     classFreqs = [6 10 12] ; %TEMP
    % updated Freqs: not inc. 6.2Hz
    classFreqs = [10 12.5 15 17.24];
seconds = 2; %2 second window
winLen = seconds*Fs; 
winFraction = 4;%125pts; %1/2 of a second %% TODO: CHANGE TO 1/4second
winShift = floor(Fs/winFraction); 
dataLimit = floor((length(fp1)-winLen)/winLen);
start = 1;
Window = cell( seconds*winFraction*dataLimit - 1, 2);
numWindows = seconds*winFraction*dataLimit;
featureH1 = cell(numWindows, length(classFreqs));
featureH2 = featureH1;
featureH3 = featureH1;
Pxx = featureH1;
PSDfeat = cell(2,1);
M1 = zeros(numWindows,2);
M2 = M1;
I1 = M1;
I2 = M1;
M3 = M1;
I3 = M1;
lx = [0 40];
ly = [1.2E-6 1.2E-6];
% % plot(lx,ly,'Color','r');
Frange = 2;
totalOps = length(classFreqs)*3*numWindows;
cont = 1;
fH = figure(1);
set(fH, 'Position', [100, 100, 1920, 1400]);
for j = 1 : length(classFreqs)
    currentFreq = classFreqs(j)
    for i = 1 : numWindows
        start = 1 + winShift*(i-1);
        winEnd = start + winLen-1;
        Window{i,1} = fp1( start : start + winLen-1 );
        Window{i,2} = fp2( start : start + winLen-1 );
        for k = 1:3 %harmonic loop:
%             figure(1);
            switch k
                case 1 % 1st harmonic
                    fp1f = ssvepFilter( Window{i,1}, Fs, k*currentFreq, Frange, 3);
                    [f, fp1fft] = get_fft_data(fp1f, Fs);
                    subplot(4,1,1);
                    [M, I] = max(fp1fft(:));
                    fp2f = ssvepFilter( Window{i,2}, Fs, k*currentFreq, Frange, 3);
                    featureH1{i,j} = [f(I),M];
                    if cont~=0
                        hold on;
                        plot(f, fp1fft),xlim(lx);
                        plot(f(I), M,'-.r*'),xlim(lx);
                        text(f(I), M, num2str(f(I)));
                        hold off;
                    end
                case 2
                    fp1f = ssvepFilter( Window{i,1}, Fs, 2*currentFreq, Frange, 3);
                    [f, fp1fft] = get_fft_data(fp1f, Fs);
                    subplot(4,1,2); 
                    [M, I] = max(fp1fft(:));
                    if cont~=0
                        hold on;
                        plot(f, fp1fft),xlim(lx);
                        plot(f(I), M,'-.r*'),xlim(lx);
                        text(f(I), M, num2str(f(I)));
                        hold off;
                    end
                    fp2f = ssvepFilter( Window{i,2}, Fs, 2*currentFreq, Frange, 3);
                    featureH2{i,j} = [f(I),M];
                case 3
                    fp1f = ssvepFilter( Window{i,1}, Fs, 2*2*currentFreq, Frange, 3);
                    [f, fp1fft] = get_fft_data(fp1f, Fs);
                    subplot(4,1,3);
                    [M, I] = max(fp1fft(:));
                    if cont~=0
                        hold on;
                        plot(f, fp1fft),xlim(lx)
                        plot(f(I), M,'-.r*'),xlim(lx);
                        text(f(I), M, num2str(f(I)));
                        hold off;
                    end
                    fp2f = ssvepFilter( Window{i,2}, Fs, 2*2*currentFreq, Frange, 3);
                    featureH3{i,j} = [f(I),M];
                otherwise
                    break;
            end
        end
        totalOps = totalOps - 3;
        if j==1 %Do only once 
%             [Pxx, ~] = pwelch( eeg_h_custom(Window{i,1},Fs,[3 30],3),[],[],Fs );
            [Pxx,ff] = pwelch( eeg_h_custom(Window{i,1},Fs,[4 30],5), 500, 250, 250, Fs);
            PSDpart1 = 9; %(<=8Hz)
            PSDend1 = 14; %(8?13Hz)
            PSDend2 = 28; %(14?27Hz)
            LPxx = 10*log10(Pxx);
            [M1(i,j), I1(i,j)] = max(LPxx(1:PSDpart1)); % 6Hz Harmonic 1
            [M2(i,j), I2(i,j)] = max(LPxx(1+PSDpart1:PSDend1)); % 6Hz Harmonic 2
            [M3(i,j), I3(i,j)] = max(LPxx(1+PSDend1:PSDend2));
            PSDfeat{j}(i,:) = [ ff(I1(i,j)), M1(i,j), ff(PSDpart1+I2(i,j)), M2(i,j), ff(PSDend1+I3(i,j)), M3(i,j)];
        end
        
        if cont~=0
            subplot(4,1,4); hold on;
            plot(ff,10*log10(Pxx));
            plot(ff(I1(i,j)), M1(i,j),'-.r*'),xlim(lx);
            text(ff(I1(i,j)), M1(i,j),num2str(ff(I1(i,j))));
            plot(ff(PSDpart1+I2(i,j)), M2(i,j),'-.r*'),xlim(lx);
            text(ff(PSDpart1+I2(i,j)), M2(i,j),num2str(ff(PSDpart1+I2(i,j))));
            plot( ff(PSDend1+I3(i,j)), M3(i,j),'-.r*'),xlim(lx);
            text( ff(PSDend1+I3(i,j)), M3(i,j),num2str(ff(PSDend1+I3(i,j))));
            hold off;
            fprintf('Operations Remaining: %d\n',totalOps);
            cont = input('continue?\n');
            clf(1);%clf(2);
        end
    end
end

%%%% Create moving window:
%%% Todo:Apply Hamming Window.
%% Organize features: (unknown value)
% testData = zeros(size(featureH1,1),32);
for j = 1:size(featureH1,1)
    testData(j,:) = [featureH1{j,1} featureH1{j,2} featureH1{j,3} featureH1{j,4}...
            featureH2{j,1} featureH2{j,2} featureH2{j,3} featureH2{j,4} ...
            featureH3{j,1} featureH3{j,2} featureH3{j,3} featureH3{j,4} ...
            PSDfeat{1}(j,:) ]; %PSDfeat{2}(j,:)
end

%% Combine Features: (Unknown):
numColns = size(freq6,2);
rows = size(freq6,1);
for i=1:rows
    freq6(i,numColns+1) = 6;
end
numColns = size(freq10,2);
rows = size(freq10,1);
for i=1:rows
    freq10(i,numColns+1) = 10;
end
numColns = size(freq12,2);
rows = size(freq12,1);
for i=1:rows
    freq12(i,numColns+1) = 12;
%%% Create moving window:
%%% Todo:Apply Hamming Window.
%% Organize features: (6Hz)
% ssFeats = zeros(numWindows, 30);
for j = 1:size(featureH1,1)
    if expectedFreq == 6
        ssFeats(j,:) = [featureH1{j,1} featureH1{j,2} featureH1{j,3} featureH1{j,4}...
            featureH2{j,1} featureH2{j,2} featureH2{j,3} featureH2{j,4} ...
            featureH3{j,1} featureH3{j,2} featureH3{j,3} featureH3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:)];
    elseif expectedFreq == 10
        ssFeats2(j,:) = [featureH1{j,1} featureH1{j,2} featureH1{j,3} featureH1{j,4}...
            featureH2{j,1} featureH2{j,2} featureH2{j,3} featureH2{j,4} ...
            featureH3{j,1} featureH3{j,2} featureH3{j,3} featureH3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:)];
    elseif expectedFreq == 12
        ssFeats3(j,:) = [featureH1{j,1} featureH1{j,2} featureH1{j,3} featureH1{j,4}...
            featureH2{j,1} featureH2{j,2} featureH2{j,3} featureH2{j,4} ...
            featureH3{j,1} featureH3{j,2} featureH3{j,3} featureH3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:)];
    elseif expectedFreq == 14
        ssFeats4(j,:) = [featureH1{j,1} featureH1{j,2} featureH1{j,3} featureH1{j,4}...
            featureH2{j,1} featureH2{j,2} featureH2{j,3} featureH2{j,4} ...
            featureH3{j,1} featureH3{j,2} featureH3{j,3} featureH3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:)];
    else
        error('wat u doin mate, you didnt put the right freq\n');
    end
end

%% Combine
numColns = size(ssFeats,2);
rows = size(ssFeats,1);
for i=1:rows
    ssFeats(i,numColns+1) = 6;
end
trainingData = [freq6; freq10; freq12];


[A, B] = trainClassifierQDA3(trainingData); %91
%% Organize all features:
% ssFeats = zeros(numWindows, 30);
for j = 1:size(featureH1,1)
    if expectedFreq == 6
        ssFeats(j,:) = [featureH1{j,1} featureH1{j,2} featureH1{j,3} featureH1{j,4}...
            featureH2{j,1} featureH2{j,2} featureH2{j,3} featureH2{j,4} ...
            featureH3{j,1} featureH3{j,2} featureH3{j,3} featureH3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:)];
    elseif expectedFreq == 10
        ssFeats2(j,:) = [featureH1{j,1} featureH1{j,2} featureH1{j,3} featureH1{j,4}...
            featureH2{j,1} featureH2{j,2} featureH2{j,3} featureH2{j,4} ...
            featureH3{j,1} featureH3{j,2} featureH3{j,3} featureH3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:)];
    elseif expectedFreq == 12
        ssFeats3(j,:) = [featureH1{j,1} featureH1{j,2} featureH1{j,3} featureH1{j,4}...
            featureH2{j,1} featureH2{j,2} featureH2{j,3} featureH2{j,4} ...
            featureH3{j,1} featureH3{j,2} featureH3{j,3} featureH3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:)];
    elseif expectedFreq == 14
        ssFeats4(j,:) = [featureH1{j,1} featureH1{j,2} featureH1{j,3} featureH1{j,4}...
            featureH2{j,1} featureH2{j,2} featureH2{j,3} featureH2{j,4} ...
            featureH3{j,1} featureH3{j,2} featureH3{j,3} featureH3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:)];
    else
        error('wat u doin mate, you didnt put the right freq\n');
    end
end

%% Combine
numColns = size(ssFeats,2);
rows = size(ssFeats,1);
for i=1:rows
    ssFeats(i,numColns+1) = 10;
end

numColns = size(ssFeats2,2);
rows = size(ssFeats2,1);
for i=1:rows
    ssFeats2(i,numColns+1) = 10;
end

numColns = size(ssFeats3,2);
rows = size(ssFeats3,1);
for i=1:rows
    ssFeats3(i,numColns+1) = 12;
end

numColns = size(ssFeats3_2,2);
rows = size(ssFeats3_2,1);
for i=1:rows
    ssFeats3_2(i,numColns+1) = 12;
end

numColns = size(ssFeats4,2);
rows = size(ssFeats4,1);
for i=1:rows
    ssFeats4(i,numColns+1) = 14;
end

SSCombine = [ssFeats;ssFeats2;ssFeats3;ssFeats3_2;ssFeats4];

%% Train Classifier
load('3classesCombinedData.mat')
[A, B] = trainClassifierQDA3(SSCombine); % 94.88%

load('4classesCombinedData.mat')
[A2, B2] = trainClassifierQDA3(SSCombine); %93.52%

% predict?
 %%% Feed me 32 features pls:
for i=1:40
    yfit(i) = A.predictFcn(testData(i,1:32));
end

for i=1:40
    yfit2(i) = A2.predictFcn(testData(i,1:32));
end

yfitcnt = (yfit == 12);
yfitcnt = sum(yfitcnt);
pct = yfitcnt / size(yfit,2);

yfitcnt2 = (yfit2 == 12);
yfitcnt2 = sum(yfitcnt2);
pct2 = yfitcnt2 / size(yfit2,2);

SSCombine = [ssFeats;ssFeats2;ssFeats3;ssFeats3_2;ssFeats4];

