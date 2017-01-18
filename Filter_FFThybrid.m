%% Import
clear all;clc;close all;
        %% All Recordings at 2KHz
% recording1 = csvread('Recording1.csv');%60Hz monitor recording
% recording1 = csvread('Recording4.csv'); %10Hz
% recording1 = csvread('Recording5.csv'); %30Hz
% recording1 = csvread('Recording7.csv'); %60Hz phone monitor
%recording1 = csvread('RecordingBaseline.csv');
% recording1 = csvread('Recording8.csv'); %LED 2x14ms delay (35.7Hz)
%recording1 = csvread('Recording9.csv'); %LED 2x7ms delay (71.43)
recording1 = csvread('31.25Hz_fadi1.csv');
%recording1 =csvread('travisphonebaseline6.csv');
%recording1= csvread('travis_baseLine_Long.csv'); % trouble reading
%baseline based on trouble reading 'numeric' field from file
%recording1=csvread('travisbaseline2.csv');
% load('hz42.mat');
%Try 6ms*2 delay
% % fp1_recording1 = BPSignal{1};
% % fp2_recording1 = BPSignal{2};
fp1_recording1 = recording1(:,1); 
fp2_recording1 = recording1(:,2);
 Fs = 250;
%Fs = 200;
seconds = 2;
start = 1;
fp1_recording1 = fp1_recording1(start*Fs:(start+seconds)*Fs,1);
clear recording1;
h=1/Fs;
L = size(fp1_recording1,1);
t = 0:h:L/Fs-h;
figure(1);
plot(t,fp1_recording1);
%% Filter (Butterworth (3.5,100Hz) +Notch at 60Hz.
figure(2);
fp1_filtered = eeg_firfilt(fp1_recording1, Fs);
plot(t, fp1_filtered);
% fp1_filtered = fp1_filtered(200:1800, 1);
%% FFT (maybe implement a smoothing filter and see how it affects measurements)

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
%% Display peak frequency:
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




