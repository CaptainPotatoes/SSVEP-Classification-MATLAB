%% Shifting Window Method
clear;clc;close all;
% load('Trial_DB1');
load('meog_t4.mat')
fp1 = Trial{1}(1:end-250,1); %ignore last second
fp2 = Trial{2}(1:end-250,1);
fpz = Trial{3}(1:end-250,1);
eyeR = Trial{4}(1:end-250,1);
% eyeL
h=1/250;
t=0:h:length(fp1)/250-h;
% markers = tD{1};
markers = trainingData{1};
Fs = SamplingRate; 
% Fs = 250;
fp1filt = eog_h_fcn(fp1,250);
fp2filt = eog_h_fcn(fp2,250);
fpzfilt = eog_h_fcn(fpz,250);
eyeRfilt = eog_h_fcn(eyeR,250);
%
figure; hold on;
plot(fp1filt,'color','r')
for i=1:length(markers)
    text(markers(i,1), fp1filt(markers(i,1)), num2str(markers(i,2)));
end
hold off;
figure; hold on;
plot(eyeRfilt,'color','c');
for i=1:length(markers)
    text(markers(i,1), eyeRfilt(markers(i,1)), num2str(markers(i,2)));
end
hold off;
figure; hold on;
plot(fpzfilt,'color','m');
for i=1:length(markers)
    text(markers(i,1), fpzfilt(markers(i,1)), num2str(markers(i,2)));
end
hold off;
figure; hold on;
plot(fp2filt,'color','y');
for i=1:length(markers)
    text(markers(i,1), fp2filt(markers(i,1)), num2str(markers(i,2)));
end
hold off;


%--ALT PLOT
close all;
f1 = figure(1);
set(f1, 'Position', [100, 100, 1600, 900]);
hold on;
plot(fp1filt,'color','r')%,ylim([-8e-4,8e-4]);
plot(eyeRfilt,'color','c');
plot(fpzfilt,'color','m');
plot(fp2filt,'color','y')%,ylim([-8e-4,8e-4]);
for i=1:length(markers)
    text(markers(i,1), fp2filt(markers(i,1)), num2str(markers(i,2)));
    if mod(i,2)==0
        text(markers(i,1), fp2filt(markers(i,1))+2.5E-5, [num2str(markers(i,1))] );
    else
        text(markers(i,1), fp2filt(markers(i,1))-2.5E-5, [num2str(markers(i,1))] );
    end
end
hold off
%% Set classes:
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
minPeakProm = 0.5E-4;
iterations = seconds*winFraction*dataLimit;
disp(iterations);
for i = 1 : iterations
    start = 1 + winShift*(i-1);
    winEnd = start + winLen-1;
    fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    Window{i,1} = fp1( start : start + winLen-1 );              % set values:
    Window{i,2} = fp2( start : start + winLen-1 );
    Window{i,3} = fpz( start : start + winLen-1 );
    Window{i,4} = eyeR( start : start + winLen-1 );
    fp1f = eogcfilt( Window{i,1} ); 
    fp2f = eogcfilt( Window{i,2} );
    fpzf = eogcfilt( Window{i,3} );
    eyeRf = eogcfilt( Window{i,4} );
    [p, l] = findpeaks(fp1f, 'MinPeakHeight',minPeakProm);
    [p1, l1] = findpeaks(fp2f, 'MinPeakHeight',minPeakProm);
    [p2, l2] = findpeaks(fpzf, 'MinPeakHeight',minPeakProm);
    hold on;
    plot(fp1f),ylim([-2.5E-4 2.5E-4]); 
    plot(fp2f),ylim([-2.5E-4 2.5E-4]);
    plot(fpzf),ylim([-2.5E-4 2.5E-4]);
    plot(eyeRf),ylim([-2.5E-4 2.5E-4]);
    plot(l,p,'-*r'); 
    plot(l1,p1,'-*c');
    plot(l2,p2,'-*y');
    hold off;
    getClass = [];
    while isempty(getClass)
        getClass = input('Enter an integer value!\n');
    end
    tY(i,1) = getClass;
    F_fp1(i,:) = featureExtractionEOG( fp1f' );
    F_fp2(i,:) = featureExtractionEOG( fp2f' );
    F_fpz(i,:) = featureExtractionEOG( fpzf' );
    F_eyeR(i,:) = featureExtractionEOG( eyeRf' );
    clf(figNum);
end
%horzcat: 32 features (8/ch)
tX = [F_fp1 F_fp2 F_fpz F_eyeR];
%


%% % % % TRAIN CLASSIFIER, KNN :: filename('knnclassification')

% Import features with extractions. 
