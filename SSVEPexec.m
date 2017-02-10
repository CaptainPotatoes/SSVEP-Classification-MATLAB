%% Clear; Load and Quick & Dirty Feature Viewing
clear all;
clc;close all;
dataRootFolder = 'C:\Users\mahmoodms\Dropbox\Public\_VCU\Yeo Lab\_SSVEP\_MATLAB-SSVEP-Classification\data';
folder{1} = ['\EOG_snap\']; folder{2} = '\SSVEP_snap\'; 
filename = 'Trial_Marc_SSVEP_31.25.mat';
% load([dataRootFolder folder{2} filename]);
% load('Marc_SSVEP_12_5Hz.mat');
load('SSVEP_8Hz_Trial1_SUBJ1.mat');
fp1 = Trial{1}; 
fp2 = Trial{2};
% Fs = 250; %Override for data w/o sampling rate.
Fs = SamplingRate
h=1/Fs;
L = size(fp1,1);
t = 0:h:L/Fs-h;
% Full Signal FFT:
    %unfilt: method 1
    figure(1)
    %{
    subplot(4,1,1);
%     fp1_fft = fft(fp1);
    fp1_fft = fft(eog_h_fcn(fp1,Fs));
    P2 = abs(fp1_fft/L);
    P1 = P2(1:(L/2+1));
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f,P1),xlim([1 75]);
    %filt: method 2
    subplot(4,1,2);
    fp1_fft_2 = fft(eeg_h_fcn3_50(fp1,Fs));
    P2 = abs(fp1_fft_2/L);
    P1 = P2(1:(L/2+1));
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f,P1),xlim([1 75]);
    %unfilt: method 1
    subplot(4,1,3);
    L2 = 2^nextpow2(L);
    fp1_fft_3 = fft(fp1, L2);
    P2 = abs(fp1_fft_3/L2);
    P1 = P2(1:(L2/2+1));
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L2/2))/L2;
    plot(f,P1),xlim([1 75]);
    %filt method 2
    subplot(4,1,4);
    fp1_fft_4 = fft(eeg_h_fcn3_50(fp1,Fs),L2);
    P2 = abs(fp1_fft_4/L2);
    P1 = P2(1:(L2/2+1));
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L2/2))/L2;
    plot(f,P1),xlim([1 75]);
    %}
fp1_f = eeg_h_fcn3_50(fp1, Fs);
% fp1_fft_2 = fft(fp1_f);
% P2 = abs(fp1_fft_2/L);
% P1 = P2(1:(L/2+1));
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
[f, P1] = get_fft_data(fp1_f, Fs);
plot(f,P1),xlim([1 50]);

n=1;
fromHz = 0;
toHz = 80;
[S{n},Fspect{n},T{n},P{n}] = spectrogram(fp1_f, 2*Fs,1*Fs,12*Fs,Fs);
figure(2)
imagesc( T{n}, Fspect{n}(Fspect{n}<toHz & Fspect{n}>fromHz), ...
             10*log10(P{n}(Fspect{n}<toHz & Fspect{n}>fromHz,:)) )
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
cb = colorbar;
ylabel(cb, 'Power (db)')
colormap(jet)
title('Channel Fp1', 'FontSize', 14)
%Pwelch:
figure(3)
hold on;
[Pxx, F] = pwelch(eeg_h_fcn3_50(fp1,Fs),[],[],250);
[Pxx2, F] = pwelch(eeg_h_fcn3_50(fp1,Fs),[],[],250);
plot(10*log10(Pxx)),xlim([0 40])
plot(10*log10(Pxx2)),xlim([0 40])
hold off;
%% Extract Features and Apply Class #
% Features: first we will use eye movement to verify (in a preceding
% window) that stimulus has changed. We can use techniques such as mean of
% all the harmonics as one of the features. We will have to put windows
% around each of our frequency hotspots. 
classes = [0 1 2 3 4];
classFreqs = [7 10 12.5 15 25];
%%% 0 = nothing/baseline noise
%   1 = first frequency (~?Hz)
%   2 = frequency (~?Hz)
%   3 = frequency (~?Hz)
%%% 4 = frequency (~?Hz).
filteredData = cell(  5,  2);
                   %Fp1, Fp2
% entire signal: @ 7Hz
figure(4); hold on;
Frange = 2;
NOrder = 3;
for i=1:length(classFreqs)
    filteredData{i,1} = ssvepFilter(fp1,Fs,classFreqs(i),Frange,NOrder);
    plot(filteredData{i,1});
    filteredData{i,2} = ssvepFilter(fp2,Fs,classFreqs(i),Frange,NOrder);
end
hold off;
classHarmonics{1} = classFreqs;
classHarmonics{2} = classFreqs*2;
classHarmonics{3} = classFreqs*3;
figure(5); hold on;
for i=1:length(classHarmonics{2})
    filteredData{length(classFreqs)+i,1} = ssvepFilter(fp1,Fs,classHarmonics{2}(i),Frange,NOrder);
    filteredData{2*length(classFreqs)+i,1} = ssvepFilter(fp1,Fs,classHarmonics{3}(i),Frange,NOrder);
    plot(filteredData{length(classFreqs)+i,1});
    filteredData{length(classFreqs)+i,2} = ssvepFilter(fp2,Fs,classHarmonics{2}(i),Frange,NOrder);
    filteredData{2*length(classFreqs)+i,1} = ssvepFilter(fp2,Fs,classHarmonics{3}(i),Frange,NOrder);
end
hold off;

%% Plot ffts:
figure(6);
fp1_P1 = cell(length(classFreqs),1);
fp1_P2 = cell(length(classFreqs),1);
fp1_fft_band = cell(3*length(classFreqs),1);
for i=1:length(classFreqs)
    subplot(length(classFreqs),1,i);
    fp1_fft_band{i} = fft(filteredData{i,1});
    L = size(filteredData{i,1},1);
    fp1_P2{i} = abs(fp1_fft_band{i}/L);
    fp1_P1{i} = fp1_P2{i}(1:(L/2+1));
    fp1_P1{i}(2:end-1) = 2*fp1_P1{i}(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f,fp1_P1{i}),xlim([1 75]);
end

figure(7);
fp1_harmonic_P1 = cell(length(classFreqs),1);
fp1_harmonic_P2 = cell(length(classFreqs),1);
% fp1_fft_band = cell(length(classFreqs),1);
for i=1:length(classFreqs)
    subplot(length(classFreqs),1,i);
    fp1_fft_band{i+length(classFreqs)} = fft(filteredData{i+length(classFreqs),1});
    L = size(filteredData{i,1},1);
    fp1_harmonic_P2{i} = abs(fp1_fft_band{i+length(classFreqs)}/L);
    fp1_harmonic_P1{i} = fp1_harmonic_P2{i}(1:(L/2+1));
    fp1_harmonic_P1{i}(2:end-1) = 2*fp1_harmonic_P1{i}(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f,fp1_harmonic_P1{i}),xlim([1 75]);
end

figure(8);
fp1_harmonic2_P1 = cell(length(classFreqs),1);
fp1_harmonic2_P2 = cell(length(classFreqs),1);
% fp1_fft_band = cell(length(classFreqs),1);
for i=1:length(classFreqs)-1
    subplot(length(classFreqs)-1,1,i);
    fp1_fft_band{i+2*length(classFreqs)} = fft(filteredData{i+2*length(classFreqs),1});
    L = size(filteredData{i,1},1);
    fp1_harmonic2_P2{i} = abs(fp1_fft_band{i+2*length(classFreqs)}/L);
    fp1_harmonic2_P1{i} = fp1_harmonic2_P2{i}(1:(L/2+1));
    fp1_harmonic2_P1{i}(2:end-1) = 2*fp1_harmonic2_P1{i}(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(f,fp1_harmonic2_P1{i}),xlim([1 75]);
end
%% Feature Extraction:
clear;close all;clc;

load('Marc_SSVEP_12_5Hz.mat');
% % Copy code from other script. 

fp1 = Trial{1}(1:end-250,1); %ignore last second
fp2 = Trial{2}(1:end-250,1);
h=1/250;
t=0:h:length(fp1)/250-h;

markers = trainingData{1};
Fs = SamplingRate; 
fp1_f = eeg_h_fcn3_50(fp1,Fs);
fp2_f = eeg_h_fcn3_50(fp2,Fs);


classFreqs = [7 10 12.5 15 25];
seconds = 2; %2 second window
winLen = seconds*Fs; 
winFraction = 2;%2.5; %1/5 of a second
winShift = floor(Fs/winFraction); 
dataLimit = floor((length(fp1)-winLen)/winLen);
start = 1;
Window = cell( seconds*winFraction*dataLimit - 1, 2);
numWindows = seconds*winFraction*dataLimit
fp1FFT = cell(length(classFreqs), numWindows);
fp2FFT = cell(length(classFreqs), numWindows);
% figure(1)

lx = [0 40];
ly = [1.2E-6 1.2E-6];
for j = 1: length(classFreqs)
    currentFreq = classFreqs(j)
    for i = 1 : numWindows
%         plot(lx,ly,'Color','r');
        start = 1 + winShift*(i-1);
        winEnd = start + winLen-1;
%         fprintf('Current index = [%d to %d]\r\n',start, winEnd);
        Window{i,1} = fp1( start : start + winLen-1 );              % set values:
        Window{i,2} = fp2( start : start + winLen-1 );
        
        fp1f = ssvepFilter( Window{i,1}, Fs, currentFreq, 3, 3);
        [fp1FFT{j,i}(:,1) , fp1FFT{j,i}(:,2)] = get_fft_data(fp1f, Fs);
        
        fp2f = ssvepFilter( Window{i,2}, Fs, currentFreq, 3, 3);
        [fp2FFT{j,i}(:,1) , fp2FFT{j,i}(:,2)] = get_fft_data(fp2f, Fs);
%         hold on;
%             plot(fp1FFT{i}(:,1), fp1FFT{i}(:,2)),xlim([0 40]);ylim([0 3E-6])
%             text(fp1FFT{i}(26,1), fp1FFT{i}(26,2),num2str(currentFreq));
%         hold off;
%         clf(1);
    end
end

% manual analysis:
manualAnalysis = input('1 to manual analyze, else quit\n');
    
if manualAnalysis == 1
   for j = 1: length(classFreqs)
       for i = 1 : numWindows
            plot(lx,ly,'Color','r');
            
            
            cont = input('1 to continue, 0 to cancel\n');
            if cont == 0
                break;
            end
       end
       if cont == 0
            break;
       end
   end
end

%%% Create moving window:
%%% Todo:Apply Hamming Window.