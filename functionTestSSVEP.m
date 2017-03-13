%%  %
    %Clear 
clear;close all;clc;
    %Import Data:
ChannelNames = {['Fp1' 'Fp2' 'Fpz' 'REye']};
% load('mssvep_t2_10_1.mat');
load('mssvep_10_1.mat');
remove = 250; % Remove final second of data.
Fs = SamplingRate;
%Import as variables and scale all to one:
remove1 = 1;
ch1 = Trial{1}(remove1:9212-remove,1);
ch2 = Trial{2}(remove1:9212-remove,1);
ch3 = Trial{3}(remove1:9212-remove,1);
% ch3 = ch3(1:2000);
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
figure; 
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
% wind = [1024 512 256 128];
[S,wfreqs] = welch_psd(ch3_f, 250, hann(1024)); 
S = S(1, :);
[S2,~] = welch_psd(ch2_f, 250, hann(1024)); 
S2 = S2(1, :);
[S3,~] = welch_psd(ch1_f, 250, hann(1024)); 
S3 = S3(1, :);

figure
hold on;
plot(wfreqs, S),xlim([1 35]);
plot(wfreqs, S2),xlim([1 35]);
plot(wfreqs, S3),xlim([1 35]);
hold off;

ch4_f = (ch1_f(1:8962)+ch2_f(1:8962)+ch3_f(1:8962))/3;
[S,wfreqs] = welch_psd(ch4_f, 250, hann(1024)); 
S = S(1, :);
figure
plot(wfreqs, S),xlim([1 35]);

%% CCA Test: (using CCA2)
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


%% 