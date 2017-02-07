%% Import
clear all;
clc;close all;
%% All Recordings at 250Hz (for now)
% load('MusaTrial1.mat')
% load('BaselineData1')
load('FadiData_31.25Hz.mat')
% load('Trial_Marc_SSVEP_15Hz.mat');
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
    subplot(4,1,1);
    fp1_fft = fft(fp1);
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
    
fp1_f = eeg_h_fcn3_50(fp1, Fs);
n=1;
fromHz = 0;
toHz = 86;
[S{n},Fspect{n},T{n},P{n}] = spectrogram(fp1_f, 2*Fs,1*Fs,12*Fs,Fs);
figure(3)
imagesc( T{n}, ...
             Fspect{n}(Fspect{n}<toHz & Fspect{n}>fromHz), ...
             10*log10(P{n}(Fspect{n}<toHz & Fspect{n}>fromHz,:)) )
%     toc
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
cb = colorbar;
ylabel(cb, 'Power (db)')
colormap(jet)
title('Channel Fp1', 'FontSize', 14)
%Pwelch:
figure(2)
[Pxx, F] = pwelch(fp1_f,[],[],250);
plot(10*log10(Pxx)),xlim([0 40])

%%Extract Features and Apply Class #
class = 
%% Create moving window:
winLen = 2*Fs; %2 seconds
winShift = Fs/5; %1/5 of a second
dataLimit = floor(length(fp1)/winLen);
start = 1;
figure(3);
Window = cell( 10*dataLimit - 1, 1);
assignedClass = cell ( 10*dataLimit - 1, 1);
hold on;
for i=1:10*dataLimit - 1
    start = 1 + 50*i;
    Window{i,1} = fp1( start : start + winLen-1 );
    assignedClass{i} = input('Enter an integer value!\n');
end
hold off;
%% Filter:

%% FFT % Implement a very specific set of windows to try and pin down signal.
%figure(3);
%fp1_fft_recording1 = fft(fp1_filtered,512);
%P2 = abs(fp1_fft_recording1/L);
%P1 = P2(1:L/2+1);
%P1(2:end-1) = 2*P1(2:end-1);
%f = Fs*(0:(L/2))/L;
%plot(f, P1), xlim([0,120])
figure(3);
f = Fs*(0:(L+10))/(L+1);
n=2^nextpow2(size(fp1_filtered,1));
fft_eeg=fft(fp1_filtered,n);
plot(f,fft_eeg)
%% Feature Extraction:



%% Classification:

% Need to make this more complex. Need to find multiple peak values
    % Looking for 35.7Hz
% [c1, c2] = find(ismember(P1, max(P1(:))))
% For 2s:
m = 2;
[b1, b2] = find(ismember(P1, max(P1(seconds*8:12*seconds))))
[c1, c2] = find(ismember(P1, max(P1(seconds*15:seconds*19))))
[d1, d2] = find(ismember(P1, max(P1(seconds*22:seconds*26))))
[e1, e2] = find(ismember(P1, max(P1(seconds*29:seconds*31))))
% FOR 1s:
% [c1, c2] = find(ismember(P1, max(P1(69:75))))
% [d1, d2] = find(ismember(P1, max(P1(76:86)))).
peak_frequency = f(1, b1)
size = P1(b1)
peak_frequency = f(1, c1)
size = P1(c1)
peak_frequency2 = f(1, d1)
size2 = P1(d1)
peak_frequency3 = f(1, e1)
size3 = P1(e1)
peakdiff = size3/size2
if abs(peak_frequency-71.43)<4%Hz
    fprintf('Success \n');
else
    fprintf('Fail \n');
end




