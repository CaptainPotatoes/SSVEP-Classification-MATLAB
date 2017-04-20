%% Shifting Window Method
clear;clc;close all;
% DATA = csvread('MarcTest4.csv');
DATA = csvread('FadiTest1.csv');
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
    cutoff(:,i) = abs(diff(customFilt(ch(:,i),Fs,[0.001 25],2)).^2);
end
f1 = figure(1); set(f1, 'Position', [100, 100, 1600, 900]); hold on; 
plot(cutoff);
for i=2:length(chf)
    if (dataTags(i) ~= dataTags(i-1)) %Signal Change
        text( i, cutoff(i,1), num2str(dataTags(i)) );
    end
end
hold off;

%% Set classes:
clear chf;
seconds = 1; %2 second window
winLen = seconds*Fs; 
winFraction = 4; %2.5; %1/4 of a second
winShift = floor( Fs/winFraction ); 
dataLimit = floor( (length(ch)-winLen)/winLen );
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
absDiffThresholdSmall = 2.75E-9;
absDiffThresholdLarge = 1.25E-8;
fprintf('Iterations: %2.2f \r\n',iterations);
chf = zeros(Fs,3);
adchf = zeros(Fs-1,3);
cnt = 1;

for i = 1 : iterations
    start = 1 + winShift*(i-1);
    winEnd = start + winLen-1;
    fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    for c = 1:3
        Window{i,c} = ch(start : start + winLen-1,c);
        chf(:,c) = customFilt(Window{i,c},Fs,[0.001 10],2);
        adchf(:,c) = abs(diff(chf(:,c))).^2;
    end
    WindowTags{i} = dataTags( start : start + winLen-1 );
    
    %%%%% TODO: CHANGE FILTER:
    
%     figure(2); hold on;
%         title('Filtered EOG Data');
%         plot(chf); hold off;
%     figure(3);hold on; 
%         plot(adchf); title('ABS(DIFF(filteredData))');hold off;
    
    thresholdCheck = find(max(adchf)>absDiffThresholdSmall);
    thresholdCheck2 = find(max(adchf)>absDiffThresholdLarge);
    %{
    if ~isempty(thresholdCheck) && isempty(thresholdCheck2)
        tagsPresent = unique(WindowTags{i},'stable')
        %%% TODO: If surpasses threshold THEN classify and store data.
%         getClass = [];
%         while isempty(getClass)
%             commandwindow;
%             getClass = input('Enter an integer value!\n');
%             if isempty(getClass)
%                 getClass = 0; 
%             end
%         end
        getClass = mode(WindowTags{i}); 
        if(getClass~=0 && getClass~=1 && getClass~=2)
            tY(cnt,1) = getClass;
            %%%% TODO: Make separate EOG Classifier/Feature Extraction. 
            F_ch1(cnt,:) = featureExtractionEOG( chf(:,1)' );
            F_ch2(cnt,:) = featureExtractionEOG( chf(:,2)' );
            F_ch3(cnt,:) = featureExtractionEOG( chf(:,3)' );
            cnt = cnt + 1;
        end
    end
    %}
    if ~isempty(thresholdCheck2)
        tagsPresent = unique(WindowTags{i},'stable')
        getClass = mode(WindowTags{i}); 
        if(getClass==1 || getClass==2)
            tY(cnt,1) = getClass;
            %%%% TODO: Make separate EOG Classifier/Feature Extraction. 
            F_ch1(cnt,:) = featureExtractionEOG( chf(:,1)' );
            F_ch2(cnt,:) = featureExtractionEOG( chf(:,2)' );
            F_ch3(cnt,:) = featureExtractionEOG( chf(:,3)' );
            cnt = cnt + 1;
        end
    end
    clf(figNum);
    clf(3);
end

%horzcat: 36 features (12/ch)
tX = [F_ch1 F_ch2 F_ch3];

clearvars -except tX tY

%% % % % TRAIN CLASSIFIER, KNN :: filename('knnclassification')
% % % % % Combine
load('fadi_eyemove_tD1.mat');
load('fadi_eyemove_tD2.mat');
load('fadi_eyemove_tD3.mat');
load('fadi_eyemove_tD4.mat');
load('fadi_eyemove_tD5.mat');
tXtY5 = [tX5 tY5];
tXtY4 = [tX4 tY4];
tXtY3 = [tX3 tY3];
tXtY2 = [tX2 tY2];
tXtY1 = [tX tY];
tXtY = [tXtY1;tXtY2;tXtY3;tXtY4;tXtY5];
clearvars -except tXtY
%{
load('fadi_blink_1.mat');
load('fadi_blink_2.mat');
load('fadi_blink_3.mat');
load('fadi_blink_4.mat');
load('fadi_blink_5.mat');

tXtY5 = [tX5 tY5];
tXtY4 = [tX4 tY4];
tXtY3 = [tX3 tY3];
tXtY2 = [tX2 tY2];
tXtY1 = [tX tY];
tXtY = [tXtY1;tXtY2;tXtY3;tXtY4;tXtY5];
%}

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