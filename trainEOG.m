%% Shifting Window Method
clear;clc;close all;
% load('Trial_DB1');
% load('meog_t4.mat')
% load('meog_t3');
load('Marc_TD_3chEOG_t1.mat');
ch1 = Trial{1}(1:end-250,1); %ignore last second
ch2 = Trial{2}(1:end-250,1);
if(length(Trial)>2) 
    ch3 = Trial{3}(1:end-250,1);
end
h=1/250;
t=0:h:length(ch1)/250-h;
markers = trainingData{1};
Fs = SamplingRate;
%% Using CSV
%{
clear;clc;close all;
data_csv = csvread('EEGSensorData_2017.04.05_11.04.57_part1.csv');
ch1 = data_csv(:,1);
ch2 = data_csv(:,2);
ch3 = data_csv(:,3);
ch4 = data_csv(:,4);
ch1f = eog_h_fcn(ch1,250);
ch2f = eog_h_fcn(ch2,250);
ch3f = eog_h_fcn(ch3,250);
ch4f = eog_h_fcn(ch4,250);
%}
ch1f = eog_h_fcn(ch1,250);
ch2f = eog_h_fcn(ch2,250);
if(length(Trial)>2) 
    ch3f = eog_h_fcn(ch3,250);
end
%%
figure; hold on;
plot(ch1f,'color','r')
for i=1:length(markers)
    text(markers(i,1), ch1f(markers(i,1)), num2str(markers(i,2)));
end
hold off;
figure; hold on;
plot(ch2f,'color','y');
for i=1:length(markers)
    text(markers(i,1), ch2f(markers(i,1)), num2str(markers(i,2)));
end
hold off;

if(length(Trial)>2) 
%     figure; hold on;plot(ch4f,'color','c');
%     for i=1:length(markers)
%         text(markers(i,1), ch4f(markers(i,1)), num2str(markers(i,2)));
%     end
%     hold off;
    figure; hold on;
    plot(ch3f,'color','m');
    for i=1:length(markers)
        text(markers(i,1), ch3f(markers(i,1)), num2str(markers(i,2)));
    end
    hold off;
end
%--ALT PLOT
close all;
f1 = figure(1);
set(f1, 'Position', [100, 100, 1600, 900]);
hold on;
plot(ch1f,'color','r')%,ylim([-8e-4,8e-4]);
plot(ch2f,'color','y')%,ylim([-8e-4,8e-4]);
if(length(Trial)>2) 
    plot(ch3f,'color','m');
end
for i=1:length(markers)
    text(markers(i,1), ch2f(markers(i,1)), num2str(markers(i,2)));
    if mod(i,2)==0
        text(markers(i,1), ch2f(markers(i,1))+2.5E-5, [num2str(markers(i,1))] );
    else
        text(markers(i,1), ch2f(markers(i,1))-2.5E-5, [num2str(markers(i,1))] );
    end
end
hold off
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
    ch1f = eogcfilt( Window{i,1} ); 
    ch2f = eogcfilt( Window{i,2} );
    if(length(Trial)>2)
        Window{i,3} = ch3( start : start + winLen-1 );
        ch3f = eogcfilt( Window{i,3} );
    end
    [p, l] = findpeaks(ch1f, 'MinPeakHeight',minPeakProm);
    [p1, l1] = findpeaks(ch2f, 'MinPeakHeight',minPeakProm);
    figure(2)
    hold on;
    plot(ch1f),ylim([-2.5E-4 2.5E-4]); 
    plot(ch2f),ylim([-2.5E-4 2.5E-4]);
    if(length(Trial)>2)
        plot(ch3f),ylim([-2.5E-4 2.5E-4]);
    end
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
    if(length(Trial)>2)
        F_ch3(i,:) = featureExtractionEOG( ch3f' );
    end
    clf(figNum);
end
%horzcat: 32 features (8/ch)
tX = [F_ch1 F_ch2 F_ch3];
%


%% % % % TRAIN CLASSIFIER, KNN :: filename('knnclassification')

% Import features with extractions. 
