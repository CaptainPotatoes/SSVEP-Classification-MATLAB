clear;clc;close all;
% LOAD TRAINING DATA: (tX, tY);
load('allEOGtD.mat');
% LOAD TEST DATA:
load('meog_t1.mat');
fp1 = Trial{1}(1:end-250,1); %ignore last second
fp2 = Trial{2}(1:end-250,1);
fpz = Trial{3}(1:end-250,1);
eyeR = Trial{4}(1:end-250,1);
Fs = SamplingRate; 
seconds = 1; %2 second window
winLen = seconds*Fs; 
winFraction = 4;%2.5; %1/4 of a second
winShift = floor(Fs/winFraction); 
dataLimit = floor((length(fp1)-winLen)/winLen);
start = 1;
numCh = 4;
Window = cell( seconds*winFraction*dataLimit - 1, numCh);
assignedClass = zeros( seconds*winFraction*dataLimit - 1, 1);
figNum = 2;
fH = figure(figNum); 
set(fH, 'Position', [100, 100, 1200, 900]);
cont = [];
for i = 1 : seconds*winFraction*dataLimit
    % TODO: Pass a 5s split window to fullHybridClassifier (similar to what
    % I did with SSVEPGUI).
    start = 1 + winShift*(i-1);
    winEnd = start + winLen-1;
    fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    Window{i,1} = fp1( start : start + winLen-1 ); % set values:
    Window{i,2} = fp2( start : start + winLen-1 );
    Window{i,3} = fpz( start : start + winLen-1 );
    Window{i,4} = eyeR( start : start + winLen-1 );
    fp1f = eog_h_fcn( Window{i,1}, Fs); 
    fp2f = eog_h_fcn( Window{i,2}, Fs);
    fpzf = eog_h_fcn( Window{i,3}, Fs);
    eyeRf = eog_h_fcn( Window{i,4}, Fs);
    hold on;
    plot(fp1f),ylim([-2.5E-4 2.5E-4]);    % Plot filtered Data. 
    plot(fp2f),ylim([-2.5E-4 2.5E-4]);
    plot(fpzf),ylim([-2.5E-4 2.5E-4]);
    plot(eyeRf),ylim([-2.5E-4 2.5E-4]);
    hold off;
    Y = fullHybridClassifier(Window{i,1}, Window{i,2}, Window{i,3}, Window{i,4}, tX, tY)
    if isempty(cont)
        cont = input('Continue? \n');
    end
    clf(2);
end

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
[S,wfreqs] = welch_estimator(ch3_f, 250, hann(1024)); 
S = S(1, :);
[S2,~] = welch_estimator(ch2_f, 250, hann(1024)); 
S2 = S2(1, :);
[S3,~] = welch_estimator(ch1_f, 250, hann(1024)); 
S3 = S3(1, :);

figure
hold on;
plot(wfreqs, S),xlim([1 35]);
plot(wfreqs, S2),xlim([1 35]);
plot(wfreqs, S3),xlim([1 35]);
hold off;

ch4_f = (ch1_f(1:8962)+ch2_f(1:8962)+ch3_f(1:8962))/3;
[S,wfreqs] = welch_estimator(ch4_f, 250, hann(1024)); 
S = S(1, :);
figure
plot(wfreqs, S),xlim([1 35]);

