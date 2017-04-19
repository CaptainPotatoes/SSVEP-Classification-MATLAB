%% Shifting Window Method
clear;clc;close all;
DATA = csvread('FadiTest5.csv');
rFB = 0; % Remove From Beginning
rFE = 0; % Remove From End
for i = 1:3
    ch(:,i) = DATA(1+rFB:end-rFE,i);
end
dataTags = DATA(1+rFB:end-rFE,4);
recordingLengthSeconds = length(ch(:,1))/250;
fprintf('Recording Length: %2.2f seconds \n',recordingLengthSeconds);
Fs = 250;
h=1/250;
t=0:h:length(ch)/250-h;
for i = 1:3
    chf(:,i) = eog_h_fcn(ch(:,i),250);
    cutoff(:,i) = abs(diff(customFilt(ch(:,i),Fs,[0.001 10],2)));
end
f1 = figure(1); set(f1, 'Position', [100, 100, 1600, 900]); hold on; 
plot(cutoff);
for i=2:length(chf)
    if (dataTags(i) ~= dataTags(i-1)) %Signal Change
        text( i, cutoff(i,1), num2str(dataTags(i)) );
    end
end
figure(2);

hold off;

%% Set classes:
seconds = 1; %2 second window
winLen = seconds*Fs; 
winFraction = 4; %2.5; %1/4 of a second
winShift = floor( Fs/winFraction ); 
dataLimit = floor( (length(ch1)-winLen)/winLen );
start = 1;
numCh = 3;
Window = cell( seconds*winFraction*dataLimit - 1, numCh );
WindowTags = cell( 1,1 ); 
assignedClass = zeros( seconds*winFraction*dataLimit - 1, 1 );
figNum = 2;
fH = figure(figNum); 
set(fH, 'Position', [100, 100, 1600, 900]);
fH2 = figure(3);
set(fH2, 'Position', [2600, 100, 1600, 900]);
YLIM = [-0.005, 0.01];
iterations = seconds*winFraction*dataLimit;
absDiffThreshold = 2.5E-5;
fprintf('Iterations: %2.2f \r\n',iterations);
cnt = 1;

for i = 1 : iterations
    start = 1 + winShift*(i-1);
    winEnd = start + winLen-1;
    fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    for c = 1:3
        Window{i,c} = ch(start : start + winLen-1,c);
    end
    WindowTags{i} = dataTags( start : start + winLen-1 );
    
    %%%%% TODO: CHANGE FILTER: 
    %%%%%%%%%% TODO: CHANGE TO FOR LOOPS.
    tch1f = eog_h_fcn(Window{i,1},250);
    tch2f = eog_h_fcn(Window{i,2},250);
    tch3f = eog_h_fcn(Window{i,3},250);
    
    dtch1f = abs(diff(tch1f));
    dtch2f = abs(diff(tch2f));
    dtch3f = abs(diff(tch3f));

    th_ch1 = abs(diff(tch1f,2));
    th_ch2 = abs(diff(tch2f,2));
    th_ch3 = abs(diff(tch3f,2));
    
    figure(2); hold on;
        plot(tch1f),ylim(YLIM); 
        plot(tch2f),ylim(YLIM); title('Filtered EOG Data');
        plot(tch3f),ylim(YLIM); hold off;
    figure(3);hold on; 
        plot(th_ch1);
        plot(th_ch2); title('ABS(DIFF6(filteredData))');
        plot(th_ch3),ylim([0 1E-4]); hold off;
    
    thresholdCheck = (max(dtch1f)>absDiffThreshold || max(dtch2f)>absDiffThreshold || max(dtch3f)>absDiffThreshold);
    if (thresholdCheck) 
        tagsPresent = unique(WindowTags{i},'stable')
        %%%% TODO: If surpasses threshold THEN classify and store data.
        getClass = [];
        while isempty(getClass)
            commandwindow;
            getClass = input('Enter an integer value!\n');
            if isempty(getClass)
                getClass = 0; 
            end
        end
        tY(cnt,1) = getClass;
        %%%% TODO: Make separate EOG Classifier/Feature Extraction. 
        F_ch1(cnt,:) = featureExtractionEOG( tch1f' );
        F_ch2(cnt,:) = featureExtractionEOG( tch2f' );
        F_ch3(cnt,:) = featureExtractionEOG( tch3f' );
        cnt = cnt + 1;
    end
    clf(figNum);
    clf(3);
end

%horzcat: 36 features (12/ch)
tX = [F_ch1 F_ch2 F_ch3];



%% % % % TRAIN CLASSIFIER, KNN :: filename('knnclassification')

% Import features with extractions. 
%         switch (dataTags(i))
%             case 0
%                 tag = 'null';
%             case 1
%                 tag = 'blink';
%             case 2
%                 tag = 'double-blink';
%             case 3
%                 tag = 'up';
%             case 4
%                 tag = 'left';
%             case 5
%                 tag = 'right';
%             case 6
%                 tag = 'down';
%             case -1
%                 tag = 'ret2center';
%         end