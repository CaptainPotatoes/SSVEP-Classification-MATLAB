%% Shifting Window Method
clear;clc;close all;
% DATA = csvread('MarcTestData.csv');
DATA = csvread('MarcTest3.csv');
rFB = 500; % Remove From Beginning
rFE = 250; % Remove From End
ch1 = DATA(1+rFB:end-rFE,1); %ignore last second
ch2 = DATA(1+rFB:end-rFE,2);
ch3 = DATA(1+rFB:end-rFE,3);
dataTags = DATA(1+rFB:end-rFE,4);
recordingLengthSeconds = length(ch1)/250;
fprintf('Recording Length: %2.2f seconds \n',recordingLengthSeconds);
Fs = 250;
h=1/250;
t=0:h:length(ch1)/250-h;
ch1f = eog_h_fcn(ch1,250);
ch2f = eog_h_fcn(ch2,250);
ch3f = eog_h_fcn(ch3,250);
f1 = figure(1); set(f1, 'Position', [100, 100, 1600, 900]); hold on;
plot(ch1f);plot(ch2f);plot(ch3f);
temp = abs(diff(ch1f));
% plot(abs(diff(ch1f)));plot(abs(diff(ch2f)));plot(abs(diff(ch3f)));
for i=2:length(ch1)
    if (dataTags(i) ~= dataTags(i-1)) %Signal Change
        switch (dataTags(i))
            case 0
                tag = 'null';
            case 1
                tag = 'blink';
            case 2
                tag = 'double-blink';
            case 3
                tag = 'up';
            case 4
                tag = 'left';
            case 5
                tag = 'right';
            case 6
                tag = 'down';
            case -1
                tag = 'ret2center';
        end
        text( i, ch1f(i), num2str(dataTags(i)));
%         text( i, temp(i), num2str(dataTags(i)));
    end
end
hold off;

%% Set classes:
seconds = 1; %2 second window
winLen = seconds*Fs; 
winFraction = 4;%2.5; %1/4 of a second
winShift = floor(Fs/winFraction); 
dataLimit = floor((length(ch1)-winLen)/winLen);
start = 1;
numCh = 4;
Window = cell( seconds*winFraction*dataLimit - 1, numCh);
assignedClass = zeros( seconds*winFraction*dataLimit - 1, 1);
figNum = 2;
fH = figure(figNum); 
set(fH, 'Position', [100, 100, 1200, 900]);
minPeakProm = 0.5E-4;
iterations = seconds*winFraction*dataLimit;
disp(iterations);
for i = 1 : iterations
    start = 1 + winShift*(i-1);
    winEnd = start + winLen-1;
    fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    Window{i,1} = ch1( start : start + winLen-1 );              % set values:
    Window{i,2} = ch2( start : start + winLen-1 );
    Window{i,3} = ch3( start : start + winLen-1 );
    %%%%% TODO: CHANGE FILTER:
    ch1f = eogcfilt( Window{i,1} ); 
    ch2f = eogcfilt( Window{i,2} );
    ch3f = eogcfilt( Window{i,3} );
    [p, l] = findpeaks(ch1f, 'MinPeakHeight',minPeakProm);
    [p1, l1] = findpeaks(ch2f, 'MinPeakHeight',minPeakProm);
    figure(2)
    hold on;
        plot(ch1f),ylim([-2.5E-4 2.5E-4]); 
        plot(ch2f),ylim([-2.5E-4 2.5E-4]);
        plot(ch3f),ylim([-2.5E-4 2.5E-4]);
        plot(l,p,'-*r'); 
        plot(l1,p1,'-*c');
    hold off;
    %%%% TODO: Make separate EOG Classifier.
    getClass = [];
    while isempty(getClass)
        commandwindow;
        getClass = input('Enter an integer value!\n');
        if isempty(getClass)
            getClass = 0; 
        end
    end
    tY(i,1) = getClass;
    F_ch1(i,:) = featureExtractionEOG( ch1f' );
    F_ch2(i,:) = featureExtractionEOG( ch2f' );
    F_ch3(i,:) = featureExtractionEOG( ch3f' );
    clf(figNum);
end
%horzcat: 32 features (8/ch)
tX = [F_ch1 F_ch2 F_ch3];
%


%% % % % TRAIN CLASSIFIER, KNN :: filename('knnclassification')

% Import features with extractions. 
