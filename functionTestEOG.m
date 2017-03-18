clear;clc;close all;
% LOAD TRAINING DATA: (tX, tY);
% load('allEOGtD.mat');
% load('tXtY_SSVEP.mat')
% LOAD TEST DATA:
% load('meog_t1.mat');
% load('mssvep_10_2.mat');
% load('mssvep_12.5_1.mat');
load('mssvep_15_1.mat');
% load('mssvep_16.6_2.mat');
% load('mssvep_t2_16_1');
fp1 = Trial{1}(1:end-250,1); %ignore last second
fp2 = Trial{2}(1:end-250,1);
fpz = Trial{3}(1:end-250,1);
eyeR = Trial{4}(1:end-250,1);
Fs = SamplingRate; 
seconds = 1; %2 second window
winLen = seconds*Fs; 
winFraction = 2;%2.5; %1/4 of a second
winShift = floor(Fs/winFraction); 
dataLimit = floor((length(fp1)-winLen)/winLen);
start = 1;
numCh = 4;
Window = cell( seconds*winFraction*dataLimit - 1, numCh);
assignedClass = zeros( seconds*winFraction*dataLimit - 1, 1);
figNum = 2;
% fH = figure(figNum); 
% set(fH, 'Position', [100, 100, 1200, 900]);
cont = [0];
Y = cell(seconds*winFraction*dataLimit,1);
% Y = zeros(seconds*winFraction*dataLimit,5);
% F = zeros(seconds*winFraction*dataLimit,53);
nS = 2;
for i = 1 : seconds*winFraction*dataLimit
    % TODO: Pass a 5s split window to fullHybridClassifier (similar to what
    % I did with SSVEPGUI).
    start = 1 + winShift*(i-1);
    winEnd = start + winLen-1;
%     fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    Window{i,1} = fp1( start : start + winLen-1 ); % set values:
    Window{i,2} = fp2( start : start + winLen-1 );
    Window{i,3} = fpz( start : start + winLen-1 );
    Window{i,4} = eyeR( start : start + winLen-1 );
    fp1f = eog_h_fcn( Window{i,1}, Fs); 
    fp2f = eog_h_fcn( Window{i,2}, Fs);
    fpzf = eog_h_fcn( Window{i,3}, Fs);
    eyeRf = eog_h_fcn( Window{i,4}, Fs);
    c1 = eegcfilt(Window{i,1});
    c2 = eegcfilt(Window{i,2});
    c3 = eegcfilt(Window{i,3});
%     hold on;
%     plot(fp1f),ylim([-2.5E-4 2.5E-4]);    % Plot filtered Data. 
%     plot(fp2f),ylim([-2.5E-4 2.5E-4]);
%     plot(fpzf),ylim([-2.5E-4 2.5E-4]);
%     plot(eyeRf),ylim([-2.5E-4 2.5E-4]);
%     hold off;
%     F(i,:) = featureExtractionSSVEP(c1,c2,c3,Fs);
    Y{i} = fullHybridClassifier(Window{i,1}, Window{i,2}, Window{i,3}, ...
        Window{i,4}, Fs, false); % boolean = EOGOnly
%     if i>nS
%         V1(((i-nS):i),:) = F(((i-nS):i),1:20);
%         V2(((i-nS):i),:) = F(((i-nS):i),21:40);
%         [A,B,r,U,V] = CCA(V1, V2);
%     end
    if isempty(cont)
        cont = input('Continue? \n');
    end
%     clf(2);
end