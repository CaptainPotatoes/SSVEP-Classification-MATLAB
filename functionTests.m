clear;close all;clc;
    %Import Data:
ChannelNames = {['Fp1' 'Fp2' 'Fpz' 'REye']};
load('mssvep_16.6_3.mat');
remove = 250; % Remove final second of data.
Fs = SamplingRate;
%Import as variables and scale all to one:

ch1 = Trial{1}(1:end-remove,1);
ch2 = Trial{2}(1:end-remove,1);
ch3 = Trial{3}(1:end-remove,1);
ch3 = ch3(1:1333);
if size(Trial,2) > 3
    ch4 = Trial{4}(1:end-remove,1);
end

flim   = [8.0 18];
winLim = [7.5 20];
N = 5;
    %Filter & Scale everything to '1'
ch1_f = scaleAbs(customFilt(ch1, Fs, flim, N));
ch2_f = scaleAbs(customFilt(ch2, Fs, flim, N));
ch3_f = scaleAbs(customFilt(ch3, Fs, flim, N));
if size(Trial,2) > 3
    ch4_f = scaleAbs(customFilt(ch4, Fs, flim, N));
end
[f,  P1]  = get_fft_data(ch1_f, Fs);
[f2, P2] = get_fft_data(ch2_f, Fs);
[f3, P3] = get_fft_data(ch3_f, Fs);
if size(Trial,2) > 3
    [f4, P4] = get_fft_data(ch4_f, Fs);
end
figure(1); 
hold on;
plot(f,  P1,'color','m'),xlim([1 35]);
plot(f2, P2,'color','c'),xlim([1 35]);
plot(f3, P3,'color','r'),xlim([1 35]);
if size(Trial,2) > 3
    plot(f4, P4,'color','b'),xlim([1 35]);
end
hold off;
title('FFT(Ch1-4)');
ylabel('|P1(f)|');
xlabel('f (Hz)');
wind = [1024 512 256 128];
[S,freqs] = welch_estimator(ch3_f, 250, hann(1024)); 
S = S(1, :);
figure
plot(freqs, S),xlim([1 35]);