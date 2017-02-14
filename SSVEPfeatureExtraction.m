%% Clear & Load Data
    %%%SSVEP FEATURE EXTRACTION
clear;close all;clc;
load('Marc_nonHair_17.24Hz_2.mat');
% expectedFreq = input('what frequency?\n');
% TEMPORARY CONT OVERRIDE:
cont = 0;
fp1 = Trial{1}(1:end-250,1); %ignore last second
fp2 = Trial{2}(1:end-250,1); 
Fs = SamplingRate;
    %%% Plot FFT (Raw)
flim = [7 30];
N = 5;
fp1_f = eeg_h_custom(fp1, Fs, flim, N);
fp2_f = eeg_h_custom(fp2, Fs, flim, N);
[f, P1] = get_fft_data(fp1_f, Fs);
[f2, P2] = get_fft_data(fp2_f, Fs);
figure; hold on;
plot(f,P1),xlim([1 35]);
plot(f2,P2),xlim([1 35]);
hold off;
Ldp = length(fp1_f); 
Lseconds = floor(length(fp1_f)/250); 
fprintf('Datapoints = %d\nSeconds = %d\n',Ldp,Lseconds);
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
    close all;
    classFreqs = [10, 12.5, 15, 17.24];
% classFreqs = 10;
seconds = 2; %2 second window
winLen = seconds*Fs; 
winFraction = 4;%125pts; %1/2 of a second %% TODO: CHANGE TO 1/4second
winShift = floor(Fs/winFraction); 
dataLimit = floor((length(fp1)-winLen)/winLen);
start = 1;
Window = cell( seconds*winFraction*dataLimit - 1, 2);
numWindows = seconds*winFraction*dataLimit;
fp1H1 = cell(numWindows, length(classFreqs));
fp1H2 = fp1H1;
fp1H3 = fp1H1;
fp2H1 = fp1H1;
fp2H2 = fp1H1;
fp2H3 = fp1H1;
Pxx1 = fp1H1;
PSDfeat = cell(size(classFreqs,2),1);
PSDfeat2 = cell(size(classFreqs,2),1);
M1 = zeros(numWindows,1);
M2 = M1;
I1 = M1;
I2 = M1;
M3 = M1;
I3 = M1;

Frange = 2;
totalOps = length(classFreqs)*3*numWindows;
% cont = 1;
fH = figure(1);
set(fH, 'Position', [100, 100, 1200, 900]);
lx = [0 40];
for j = 1 : length(classFreqs)
    currentFreq = classFreqs(j)
    for i = 1 : numWindows
        start = 1 + winShift*(i-1);
        winEnd = start + winLen-1;
        Window{i,1} = fp1( start : start + winLen-1 );
        Window{i,2} = fp2( start : start + winLen-1 );
        for k = 1:3 %harmonic loop:
            switch k
                case 1 % 1st harmonic
                    % CH1
                    fp1f = ssvepFilter( Window{i,1}, Fs, (2^(k-1))*currentFreq, Frange, 3);
                    [f, fp1fft] = get_fft_data(fp1f, Fs);
                    [M, I] = max(fp1fft(:));
                    % CH2
                    fp2f = ssvepFilter( Window{i,2}, Fs, (2^(k-1))*currentFreq, Frange, 3);
                    [f2, fp2fft] = get_fft_data(fp2f, Fs); 
                    [M2, I2] = max(fp2fft(:));
                    %Plot Both: (SUBPLOT 1 of 4)
                    subplot(4,1,1);
                    if cont~=0
                        hold on;
                        plot(f, fp1fft),xlim(lx);
                        plot(f2, fp2fft),xlim(lx);
                        plot(f(I), M,'-.r*'),xlim(lx);
                        plot(f2(I2), M2,'-.b^'),xlim(lx);
                        text(f(I), M, num2str(f(I)));
                        text(f2(I2), M2, num2str(f2(I2)));
                        hold off;
                    end
                    fp1H1{i,j} = [f(I),M];
                    fp2H1{i,j} = [f2(I2),M2];
                case 2
                    % CH1
                    fp1f = ssvepFilter( Window{i,1}, Fs, (2^(k-1))*currentFreq, Frange, 3);
                    [f, fp1fft] = get_fft_data(fp1f, Fs);
                    [M, I] = max(fp1fft(:));
                    % CH2
                    fp2f = ssvepFilter( Window{i,2}, Fs, (2^(k-1))*currentFreq, Frange, 3);
                    [f2, fp2fft] = get_fft_data(fp2f, Fs); 
                    [M2, I2] = max(fp2fft(:));
                    %Plot Both:
                    subplot(4,1,2);
                    if cont~=0
                        hold on;
                        plot(f, fp1fft),xlim(lx);
                        plot(f2, fp2fft),xlim(lx);
                        plot(f(I), M,'-.r*'),xlim(lx);
                        plot(f2(I2), M2,'-.b^'),xlim(lx);
                        text(f(I), M, num2str(f(I)));
                        text(f2(I2), M2, num2str(f2(I2)));
                        hold off;
                    end
                    fp1H2{i,j} = [f(I),M];
                    fp2H2{i,j} = [f2(I2),M2];
                case 3
                    % CH1
                    fp1f = ssvepFilter( Window{i,1}, Fs, (2^(k-1))*currentFreq, Frange, 3);
                    [f, fp1fft] = get_fft_data(fp1f, Fs);
                    [M, I] = max(fp1fft(:));
                    % CH2
                    fp2f = ssvepFilter( Window{i,2}, Fs, (2^(k-1))*currentFreq, Frange, 3);
                    [f2, fp2fft] = get_fft_data(fp2f, Fs); 
                    [M2, I2] = max(fp2fft(:));
                    %Plot Both:
                    subplot(4,1,3);
                    if cont~=0
                        hold on;
                        plot(f, fp1fft),xlim(lx);
                        plot(f2, fp2fft),xlim(lx);
                        plot(f(I), M,'-.r*'),xlim(lx);
                        plot(f2(I2), M2,'-.b^'),xlim(lx);
                        text(f(I), M, num2str(f(I)));
                        text(f2(I2), M2, num2str(f2(I2)));
                        hold off;
                    end
                    fp1H3{i,j} = [f(I),M];
                    fp2H3{i,j} = [f2(I2),M2];
                otherwise
                    break;
            end
        end
        totalOps = totalOps - 3;
%         [Pxx0,Psf] = pwelch( eeg_h_custom(Window{i,1}, Fs, [5 30], 5), 500, 250, 250, Fs );
%         [Pxx02,Psf2] = pwelch( eeg_h_custom(Window{i,2}, Fs, [5 30], 5), 500, 250, 250, Fs );
%         LPxx0 = 10*log10(Pxx0);
%         LPxx02 = 10*log10(Pxx02);
        switch j
            case 1 %10Hz
                 %CH1
                [Pxx1,ff] = pwelch( ssvepFilter(Window{i,1},Fs, currentFreq, Frange, 3), 500, 250, 250, Fs );
                LPxx1 = 10*log10(Pxx1);
                [M1(i), I1(i)] = max(LPxx1(9:13)); %[from 8Hz?12Hz]
                PSDfeat{j}(i,:) = [ ff(9+I1(i)), M1(i) ];
                 %CH2:
                [Pxx2,ff2] = pwelch( ssvepFilter(Window{i,2},Fs, currentFreq, Frange, 3), 500, 250, 250, Fs ); 
                LPxx2 = 10*log10(Pxx2);
                [PM2, PI2] = max(LPxx2(9:13));
                PSDfeat2{j}(i,:) = [ ff2(9+PI2), PM2 ];
                %Plot Both
                if cont~=0
                    subplot(4,1,4); hold on;
%                     plot(Psf,LPxx0);%plot full current spectrum
                    
                    plot(ff, LPxx1);%CH1
                    plot(ff(8+I1(i)), M1(i),'-.b*'),xlim(lx);
                    text(ff(8+I1(i)), M1(i),num2str(ff(8+I1(i))));
                    
                    plot(ff2,LPxx2);%CH2
                    plot(ff2(8+PI2), PM2,'-.g*'),xlim(lx);
                    text(ff2(8+PI2), PM2,num2str(ff2(8+PI2)));
                    hold off;
                end
            case 2 %12.5Hz
                 %CH1
                [Pxx1,ff] = pwelch( ssvepFilter(Window{i,1},Fs, currentFreq, Frange, 3), 500, 250, 250, Fs );
                LPxx1 = 10*log10(Pxx1);
                [M1(i), I1(i)] = max(LPxx1(12:14)); %[from 11Hz?13Hz]
                PSDfeat{j}(i,:) = [ ff(12+I1(i)), M1(i) ];
                 %CH2:
                [Pxx2,ff2] = pwelch( ssvepFilter(Window{i,2},Fs, currentFreq, Frange, 3), 500, 250, 250, Fs ); 
                LPxx2 = 10*log10(Pxx2);
                [PM2, PI2] = max(LPxx2(12:14));
                PSDfeat2{j}(i,:) = [ ff2(12+PI2), PM2 ];
                %Plot Both
                if cont~=0
                    subplot(4,1,4); hold on;
%                     plot(Psf,LPxx0);%plot full current spectrum
                    
                    plot(ff, LPxx1); %CH1
                    plot(ff(11+I1(i)), M1(i),'-.b*'),xlim(lx);
                    text(ff(11+I1(i)), M1(i),num2str(ff(11+I1(i))));
                    
                    plot(ff2,LPxx2);%CH2
                    plot(ff2(11+PI2), PM2,'-.g*'),xlim(lx);
                    text(ff2(11+PI2), PM2,num2str(ff2(11+PI2)));
                    
                    hold off;
                end
            case 3 %15
                 %CH1
                [Pxx1,ff] = pwelch( ssvepFilter(Window{i,1},Fs, currentFreq, Frange, 3), 500, 250, 250, Fs );
                LPxx1 = 10*log10(Pxx1);
                [M1(i), I1(i)] = max(LPxx1(15:17)); %[from 14Hz?16Hz]
                PSDfeat{j}(i,:) = [ ff(15+I1(i)), M1(i) ];
                 %CH2:
                [Pxx2,ff2] = pwelch( ssvepFilter(Window{i,2},Fs, currentFreq, Frange, 3), 500, 250, 250, Fs ); 
                LPxx2 = 10*log10(Pxx2);
                [PM2, PI2] = max(LPxx2(15:17));
                PSDfeat2{j}(i,:) = [ ff2(15+PI2), PM2 ];
                if cont~=0
                    subplot(4,1,4); hold on;
%                     plot(Psf,LPxx0);%plot full current spectrum
                    
                    plot(ff, LPxx1); %CH1
                    plot(ff(14+I1(i)), M1(i),'-.b*'),xlim(lx);
                    text(ff(14+I1(i)), M1(i),num2str(ff(14+I1(i))));
                    
                    plot(ff2,LPxx2);%CH2
                    plot(ff2(14+PI2), PM2,'-.g*'),xlim(lx);
                    text(ff2(14+PI2), PM2,num2str(ff2(14+PI2)));
                    hold off;
                end
            case 4 %17.24
                 %CH1
                [Pxx1,ff] = pwelch( ssvepFilter(Window{i,1},Fs, currentFreq, Frange, 3), 500, 250, 250, Fs );
                LPxx1 = 10*log10(Pxx1);
                [M1(i), I1(i)] = max(LPxx1(18:28)); %[from 17Hz?27Hz]
                PSDfeat{j}(i,:) = [ ff(18+I1(i)), M1(i) ];
                 %CH2:
                [Pxx2,ff2] = pwelch( ssvepFilter(Window{i,2},Fs, currentFreq, Frange, 3), 500, 250, 250, Fs ); 
                LPxx2 = 10*log10(Pxx2);
                [PM2, PI2] = max(LPxx2(18:28));
                PSDfeat2{j}(i,:) = [ ff2(18+PI2), PM2 ];
                if cont~=0
                    subplot(4,1,4); hold on;
%                     plot(Psf,LPxx0);%plot full current spectrum
                    
                    plot(ff, LPxx1); %CH1
                    plot(ff(17+I1(i)), M1(i),'-.b*'),xlim(lx);
                    text(ff(17+I1(i)), M1(i),num2str(ff(17+I1(i))));
                    
                    plot(ff2,LPxx2);%CH2
                    plot(ff2(17+PI2), PM2,'-.g*'),xlim(lx);
                    text(ff2(17+PI2), PM2,num2str(ff2(17+PI2)));
                    hold off;
                end
            otherwise
                break;
        end
        if cont~=0
            fprintf('Operations Remaining: %d\n',totalOps);
            cont = input('continue?\n');
            clf(1);
        end
        %{
        if j==1 %Do only once 
%             [Pxx,ff] = pwelch( eeg_h_custom(Window{i,1},Fs,[4 30],5), 500, 250, 250, Fs );
            [Pxx1,ff] = pwelch( ssvepFilter(Window{i,1},Fs, currentFreq, Frange, 3), 500, 250, 250, Fs );
            PSDstart = 9;
            PSDpart1 = 14; %(<=8Hz)
            PSDend1 = 17; %(8?13Hz)
            PSDend2 = 28; %(14?27Hz)
            LPxx = 10*log10(Pxx1);
            [M1(i,j), I1(i,j)] = max(LPxx(PSDstart:PSDpart1)); % [
            [M2(i,j), I2(i,j)] = max(LPxx(1+PSDpart1:PSDend1)); % 6Hz Harmonic 2
            [M3(i,j), I3(i,j)] = max(LPxx(1+PSDend1:PSDend2));
            PSDfeat(i,:) = [ ff(I1(i,j)), M1(i,j), ff(PSDpart1+I2(i,j)), M2(i,j), ff(PSDend1+I3(i,j)), M3(i,j)];
        end
        if cont~=0
            subplot(4,1,4); hold on;
            plot(ff,10*log10(Pxx1));
            plot(ff(PSDstart+I1(i,j)), M1(i,j),'-.b*'),xlim(lx);
            text(ff(PSDstart+I1(i,j)), M1(i,j),num2str(ff(PSDstart+I1(i,j))));
            plot(ff(PSDpart1+I2(i,j)), M2(i,j),'-.b*'),xlim(lx);
            text(ff(PSDpart1+I2(i,j)), M2(i,j),num2str(ff(PSDpart1+I2(i,j))));
            plot( ff(PSDend1+I3(i,j)), M3(i,j),'-.b*'),xlim(lx);
            text( ff(PSDend1+I3(i,j)), M3(i,j),num2str(ff(PSDend1+I3(i,j))));
            hold off;
            fprintf('Operations Remaining: %d\n',totalOps);
            cont = input('continue?\n');
            clf(1);%clf(2);
        end
        %}
    end
end

%%%% Create moving window:
%%% Todo:Apply Hamming Window.
%% Organize features: (unknown value)
clearvars -except fp1H1 fp1H2 fp1H3 fp2H1 fp2H2 fp2H3 PSDfeat PSDfeat2
testData = zeros(size(fp1H1,1),64);
for j = 1:size(fp1H1,1)
    testData(j,:) = [fp1H1{j,1} fp1H1{j,2} fp1H1{j,3} fp1H1{j,4}...
            fp1H2{j,1} fp1H2{j,2} fp1H2{j,3} fp1H2{j,4} ...
            fp1H3{j,1} fp1H3{j,2} fp1H3{j,3} fp1H3{j,4} ...
            fp2H1{j,1} fp2H1{j,2} fp2H1{j,3} fp2H1{j,4}...
            fp2H2{j,1} fp2H2{j,2} fp2H2{j,3} fp2H2{j,4} ...
            fp2H3{j,1} fp2H3{j,2} fp2H3{j,3} fp2H3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:) PSDfeat{3}(j,:) PSDfeat{4}(j,:) ...
            PSDfeat2{1}(j,:) PSDfeat2{2}(j,:) PSDfeat2{3}(j,:) PSDfeat2{4}(j,:)];
end
clearvars -except testData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combine same freqs: (HOWEVER MANY THERE ARE!)
    freq10 = [f10_2;f10_3];
    freq12 = [f12_3;f12_4];
    freq15 = [f15_1;f15_2;f15_3];
    freq17 = [f17_1;f17_2;f17_3;f17_4;f17_5];
    
%% Combine ALL Data:
    freq10 = [f10_1;f10_2;f10_3;f10_4;f10_5;f10_6;f10_7];
    freq12 = [f12_1;f12_2;f12_3;f12_4;f12_5];
    freq15 = [f15_1;f15_2;f15_3;f15_4;f15_5];
    freq17 = [f17_1;f17_2;f17_3;f17_4;f17_5];
%% Label and Combine all Training Data
numColns = size(freq10,2);
rows = size(freq10,1);
for i=1:rows
    freq10(i,numColns+1) = 10;
end

numColns = size(freq12,2);
rows = size(freq12,1);
for i=1:rows
    freq12(i,numColns+1) = 12;
end

numColns = size(freq15,2);
rows = size(freq15,1);
for i=1:rows
    freq15(i,numColns+1) = 15;
end

numColns = size(freq17,2);
rows = size(freq17,1);
for i=1:rows
    freq17(i,numColns+1) = 17;
end

combinedTrainingData = [freq10;freq12;freq15;freq17];
% cALLTD = [freq10;freq12;freq15;freq17];
%%% Create moving window:
%%% Todo:Apply Hamming Window.

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
for j = 1:size(fp1H1,1)
    if expectedFreq == 6
        ssFeats(j,:) = [fp1H1{j,1} fp1H1{j,2} fp1H1{j,3} fp1H1{j,4}...
            fp1H2{j,1} fp1H2{j,2} fp1H2{j,3} fp1H2{j,4} ...
            fp1H3{j,1} fp1H3{j,2} fp1H3{j,3} fp1H3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:)];
    elseif expectedFreq == 10
        ssFeats2(j,:) = [fp1H1{j,1} fp1H1{j,2} fp1H1{j,3} fp1H1{j,4}...
            fp1H2{j,1} fp1H2{j,2} fp1H2{j,3} fp1H2{j,4} ...
            fp1H3{j,1} fp1H3{j,2} fp1H3{j,3} fp1H3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:)];
    elseif expectedFreq == 12
        ssFeats3(j,:) = [fp1H1{j,1} fp1H1{j,2} fp1H1{j,3} fp1H1{j,4}...
            fp1H2{j,1} fp1H2{j,2} fp1H2{j,3} fp1H2{j,4} ...
            fp1H3{j,1} fp1H3{j,2} fp1H3{j,3} fp1H3{j,4} ...
            PSDfeat{1}(j,:) PSDfeat{2}(j,:)];
    elseif expectedFreq == 14
        ssFeats4(j,:) = [fp1H1{j,1} fp1H1{j,2} fp1H1{j,3} fp1H1{j,4}...
            fp1H2{j,1} fp1H2{j,2} fp1H2{j,3} fp1H2{j,4} ...
            fp1H3{j,1} fp1H3{j,2} fp1H3{j,3} fp1H3{j,4} ...
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

