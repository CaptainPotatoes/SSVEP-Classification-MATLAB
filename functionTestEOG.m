clear;clc;close all;
% LOAD TRAINING DATA: (tX, tY);
% load('allEOGtD.mat');
% load('tXtY_SSVEP.mat')
% LOAD TEST DATA:
% load('meog_t1.mat');
% load('mssvep_10_2.mat');
% load('mssvep_12.5_1.mat');
% load('mssvep_15_1.mat');
load('mssvep_16.6_3.mat');
% load('mssvep_10_5');
% load('Marc_TEST_10.mat');
% load('meog_t3')
removeStart=500;
ch1 = Trial{1}(1+removeStart:end-250,1); %ignore last second
ch2 = Trial{2}(1+removeStart:end-250,1);
ch3 = Trial{3}(1+removeStart:end-250,1);
ch4 = Trial{4}(1+removeStart:end-250,1);
Fs = SamplingRate; 
range = 250:60:2500;
% range = 500:60:2000;
Window = cell(size(range,2),4);
Y = cell(size(range,2),1);
cont = [];
% Y = zeros(seconds*winFraction*dataLimit,5);
for i = 1:size(range,2)
    start = 1;
    winEnd = start + (range(i)-1);
    fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    fprintf('Seconds Elapsed = [%1.2f]\r\n',winEnd/250);
    Window{i,1} = ch1( start : winEnd ); % set values:
    Window{i,2} = ch2( start : winEnd );
    Window{i,3} = ch3( start : winEnd );
    Window{i,4} = ch4( start : winEnd );
    [Y{i},F{i}] = fHC(Window{i,1}, Window{i,2}, Window{i,3}, ...
        Window{i,4}, Fs, false);
    O10 = F{i}(1:4,:)
    O12 = F{i}(5:8,:)
    O15 = F{i}(9:12,:)
    O16 = F{i}(13:16,:)
    B10 = sum(O10(:)==0)
    B12 = sum(O12(:)==0)
    B15 = sum(O15(:)==0)
    B16 = sum(O16(:)==0)
    if isempty(cont)
        commandwindow;
        cont = input('Continue? \n');
    end
end
%{
seconds = 2; %2 second window
winLen = seconds*Fs; 
winFraction = 2;%2.5; %1/4 of a second
winShift = floor(Fs/winFraction); 
dataLimit = floor((length(ch1)-winLen)/winLen);
start = 1;
numCh = 4;
Window = cell( seconds*winFraction*dataLimit - 1, numCh);
assignedClass = zeros( seconds*winFraction*dataLimit - 1, 1);
figNum = 2;
% fH = figure(figNum); 
% set(fH, 'Position', [100, 100, 1200, 900]);
cont = [];
Y = cell(seconds*winFraction*dataLimit,1);
for i = 1 : seconds*winFraction*dataLimit
    % TODO: Pass a 5s split window to fullHybridClassifier (similar to what
    % I did with SSVEPGUI).
    start = 1 + winShift*(i-1);
    winEnd = start + winLen-1;
%     fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    Window{i,1} = ch1( start : start + winLen-1 ); % set values:
    Window{i,2} = ch2( start : start + winLen-1 );
    Window{i,3} = ch3( start : start + winLen-1 );
    Window{i,4} = ch4( start : start + winLen-1 );
    ch1f = eog_h_fcn( Window{i,1}, Fs); 
    ch2f = eog_h_fcn( Window{i,2}, Fs);
    ch3f = eog_h_fcn( Window{i,3}, Fs);
    ch4f = eog_h_fcn( Window{i,4}, Fs);
    c1 = eegcfilt(Window{i,1});
    c2 = eegcfilt(Window{i,2});
    c3 = eegcfilt(Window{i,3});
%     hold on;
%         plot(ch1f),ylim([-2.5E-4 2.5E-4]);    % Plot filtered Data. 
%         plot(ch2f),ylim([-2.5E-4 2.5E-4]);
%         plot(ch3f),ylim([-2.5E-4 2.5E-4]);
%         plot(ch4f),ylim([-2.5E-4 2.5E-4]);
%     hold off;
%     F(i,:) = featureExtractionSSVEP(c1,c2,c3,Fs);
    Y{i} = fullHybridClassifier(Window{i,1}, Window{i,2}, Window{i,3}, ...
        Window{i,4}, Fs, false); % boolean = EOGOnly
    OUTPUT = Y{i}'
    if isempty(cont)
        commandwindow;
        cont = input('Continue? \n');
    end
end
%}

% % EOF


