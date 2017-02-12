%% Clear & Load Data
clear;close all;clc;
%% 
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
filename = 'Dad_X1_12Hz.mat';
    disp(filename);
load([dataRootFolder folder{3} filename]);
% load('Marc_SSVEP_12_5Hz.mat');
fp1 = Trial{1}(1:end-250,1); %ignore last second
fp2 = Trial{2}(1:end-250,1);
Fs = SamplingRate;
%%% Plot FFT (Raw)
f = [3 40];
N = 3;
fp1_f = eeg_h_custom(fp1, Fs, f, N);
fp2_f = eeg_h_custom(fp2, Fs, f, N);
[f, P1] = get_fft_data(fp1_f, Fs);
[f2, P2] = get_fft_data(fp2_f, Fs);
figure; hold on;
plot(f,P1),xlim([1 50]);
plot(f2,P2),xlim([1 50]);
hold off;
Lseconds = length(fp1_f);
h=1/250;
t=0:h:length(fp1)/250-h;
markers = trainingData{1};
%Spect:
figure;
[~, Fspect, T, P] = spectrogram(fp1_f(500:end-500), 5*Fs,4*Fs,10*Fs,Fs);
imagesc(T, Fspect(Fspect<50 & Fspect>1), 10*log10(P(Fspect<50 & Fspect>1,:)));
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
cb = colorbar;
ylabel(cb, 'Power (db)')
colormap(jet)
title('Channel Fp1', 'FontSize', 14)
figure;
[~, Fspect, T, P] = spectrogram(fp2_f, 5*Fs,4*Fs,10*Fs,Fs);
imagesc(T, Fspect(Fspect<50 & Fspect>1), 10*log10(P(Fspect<50 & Fspect>1,:)));
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
cb = colorbar;
ylabel(cb, 'Power (db)')
colormap(jet)
title('Channel Fp1', 'FontSize', 14)
%%%Analysis
    close all;
% classFreqs = [6 10 13? 15??];
classFreqs = [6 10];
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
% PSDfeat = zeros(numWindows,4);
M1 = zeros(numWindows,2);
M2 = M1;
I1 = M1;
I2 = M1;
% figure(1)
lx = [0 40];
ly = [1.2E-6 1.2E-6];
% % plot(lx,ly,'Color','r');
Frange = 2;
totalOps = length(classFreqs)*3*numWindows
cont = 1;
for j = 1 : length(classFreqs)
    currentFreq = classFreqs(j)
    for i = 1 : numWindows
        start = 1 + winShift*(i-1);
        winEnd = start + winLen-1;
        Window{i,1} = fp1( start : start + winLen-1 );
        Window{i,2} = fp2( start : start + winLen-1 );
        for k = 1:3 %harmonic loop:
            figure(1);
            switch k
                case 1 % 1st harmonic
                    fp1f = ssvepFilter( Window{i,1}, Fs, k*currentFreq, Frange, 3);
                    [f, fp1fft] = get_fft_data(fp1f, Fs);
                    subplot(3,1,1);
                    [M, I] = max(fp1fft(:));
                    fp2f = ssvepFilter( Window{i,2}, Fs, k*currentFreq, Frange, 3);
                    featureH1{i,j} = [f(I),M];
                    if cont~=0
                        hold on;
                        plot(f, fp1fft),xlim(lx);
                        plot(f(I), M,'-.r*'),xlim(lx);
                        hold off;
                    end
                case 2
                    fp1f = ssvepFilter( Window{i,1}, Fs, 2*currentFreq, Frange, 3);
                    [f, fp1fft] = get_fft_data(fp1f, Fs);
                    subplot(3,1,2); 
                    [M, I] = max(fp1fft(:));
                    if cont~=0
                        hold on;
                        plot(f, fp1fft),xlim(lx);
                        plot(f(I), M,'-.r*'),xlim(lx);
                        hold off;
                    end
                    fp2f = ssvepFilter( Window{i,2}, Fs, 2*currentFreq, Frange, 3);
                    featureH2{i,j} = [f(I),M];
                case 3
                    fp1f = ssvepFilter( Window{i,1}, Fs, 2*2*currentFreq, Frange, 3);
                    [f, fp1fft] = get_fft_data(fp1f, Fs);
                    subplot(3,1,3);
                    [M, I] = max(fp1fft(:));
                    if cont~=0
                        hold on;
                        plot(f, fp1fft),xlim(lx)
                        plot(f(I), M,'-.r*'),xlim(lx);
                        hold off;
                    end
                    fp2f = ssvepFilter( Window{i,2}, Fs, 2*2*currentFreq, Frange, 3);
                    featureH3{i,j} = [f(I),M];
                otherwise
                    break;
            end
        end
        % TODO: PSD for each window
        
        totalOps = totalOps - 3;
        if j==1 %Do only once 
            [Pxx, ~] = pwelch( eeg_h_custom(Window{i,1},Fs,[3 30],3),[],[],Fs );
            fpsd = 1:length(Pxx);
            PSDsplit = 11;
            [M1(i,j), I1(i,j)] = max(Pxx(1:PSDsplit));
            [M2(i,j), I2(i,j)] = max(Pxx(1+PSDsplit:end));
            PSDfeat{j}(i,:) = [ fpsd(I1(i,j)), 10*log10(M1(i,j)), fpsd(PSDsplit+I2(i,j)), 10*log10(M2(i,j)) ];
        end
        if j==2
            [Pxx, ~] = pwelch( eeg_h_custom(Window{i,1},Fs,[3 30],3),[],[],Fs );
            fpsd = 1:length(Pxx);
            PSDsplit = 14;
            [M1(i,j), I1(i,j)] = max(Pxx(1:PSDsplit));
            [M2(i,j), I2(i,j)] = max(Pxx(1+PSDsplit:end));
            PSDfeat{j}(i,:) = [ fpsd(I1(i,j)), 10*log10(M1(i,j)), fpsd(PSDsplit+I2(i,j)), 10*log10(M2(i,j)) ];
        end
        if cont~=0
            figure(2); hold on;
            plot(10*log10(Pxx));
            plot(fpsd(I1(i,j)), 10*log10(M1(i,j)),'-.r*'),xlim(lx);
            text(fpsd(I1(i,j)), 10*log10(M1(i,j)),num2str(fpsd(I1(i,j))));
            plot(fpsd(PSDsplit+I2(i,j)), 10*log10(M2(i,j)),'-.r*'),xlim(lx);
            text(fpsd(PSDsplit+I2(i,j)), 10*log10(M2(i,j)),num2str(fpsd(PSDsplit+I2(i,j))));
            hold off;
            fprintf('Operations Remaining: %d\n',totalOps);
            cont = input('continue?\n');
            clf(1);clf(2);
        end
    end
end

%%% Create moving window:
%%% Todo:Apply Hamming Window.
%% Organize features: (6Hz)
for j = 1:size(featureH1,1)
    ssFeats(j,:)  = [featureH1{j,1} featureH1{j,2} featureH2{j,1} featureH2{j,2} featureH3{j,1} featureH3{j,2} PSDfeat{1}(j,:) PSDfeat{2}(j,:)]; %
end

%% Organize features: (10Hz)
for j = 1:size(featureH1,1)
    ssFeats2(j,:) = [featureH1{j,1} featureH1{j,2} featureH2{j,1} featureH2{j,2} featureH3{j,1} featureH3{j,2} PSDfeat{1}(j,:) PSDfeat{2}(j,:)]; 
end

%% Organize features: (??Hz)
for j = 1:size(featureH1,1)
    ssFeats3(j,:) = [featureH1{j,1} featureH1{j,2} featureH2{j,1} featureH2{j,2} featureH3{j,1} featureH3{j,2} PSDfeat{1}(j,:) PSDfeat{2}(j,:)]; 
end


%% Combine
 numColns = size(ssFeats,2);
 rows = size(ssFeats,1);
for i=1:rows
    ssFeats(i,numColns+1) = 6;
end

 numColns = size(ssFeats2,2);
 rows = size(ssFeats2,1);
for i=1:rows
    ssFeats2(i,numColns+1) = 10;
end

SSCombine = [ssFeats;ssFeats2];
