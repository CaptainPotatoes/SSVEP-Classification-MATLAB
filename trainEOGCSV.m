%% Shifting Window Method
clear;clc;close all;
% DATA = csvread('MarcTest0.csv');
DATA = csvread('MarcDiff2ch_T1.csv');
% DATA2 = csvread('MarcDiff2ch_T2.csv');
% DATA3 = csvread('MarcDiff2ch_T3.csv');
% DATA4 = csvread('MarcDiff2ch_T4.csv');
% DATA5 = csvread('MarcDiff2ch_T5.csv');
% DATA = [DATA1;DATA2;DATA3;DATA4;DATA5];
rFB = 0; % Remove From Beginning
rFE = 0; % Remove From End
numCh = 2;
for i = 1:numCh
    ch(:,i) = DATA(1+rFB:end-rFE,i);
end
THRESHOLD1 = 2.85E-4;
THRESHOLD2 = 3.6E-3;
THRESHOLD3 = -0.5E-3;
dataTags = DATA(1+rFB:end-rFE,3);
recordingLengthSeconds = length(ch(:,1))/250;
fprintf('Recording Length: %2.2f seconds \n',recordingLengthSeconds);
Fs = 250;
h=1/250;
for i = 1:numCh
    chf(:,i) = customFilt(ch(:,i),Fs,[1 10],2);
%     cutoff(:,i) = abs(diff(customFilt(ch(:,i),Fs,[0.2 15],2)).^2);
    cutoff(:,i) = diff(chf(:,i));
end
f1 = figure(1); set(f1, 'Position', [100, 100, 1600, 900]); hold on;
title('EOG Signal Amplitude');
xlabel('Time (s)'); 
ylabel('EOG Signal Amplitude (mV)');
s = length(ch)/250;
h=1/250;
t=0:h:(s-h);
plot(t,chf),xlim([0,max(t)]);hold on;legend('EOG Channel 1','EOG Channel 2')
rl1 = refline(0,THRESHOLD1); rl1.Color = 'r';
rl2 = refline(0,THRESHOLD2); rl2.Color = 'm';
rl3 = refline(0,THRESHOLD3); rl3.Color = 'c';
for i=2:length(chf)
    if (dataTags(i) ~= dataTags(i-1)) %Signal Change
        text( t(i), chf(i,1), num2str(dataTags(i)) );
    end
end
f2 = figure(2); set(f2, 'Position', [2600, 100, 1600, 900]); hold on; 
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
Window = cell( seconds*winFraction*dataLimit - 1, numCh );
WindowTags = cell( 1,1 ); 
assignedClass = zeros( seconds*winFraction*dataLimit - 1, 1 );
figNum = 2;
fH = figure(figNum); 
set(fH, 'Position', [100, 100, 1200, 600]);
fH2 = figure(3);
set(fH2, 'Position', [1940, 100, 900, 600]);
YLIM = [-0.005, 0.01];
iterations = seconds*winFraction*dataLimit;
fprintf('Iterations: %2.2f \r\n',iterations);
chf = zeros(Fs,3);
adchf = zeros(Fs-1,3);
cnt = 1;

for i = 1 : iterations
    start = 1 + winShift*(i-1);
    winEnd = start + winLen-1;
    fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    for c = 1:numCh
        Window{i,c} = ch(start : start + winLen-1,c);
        chf(:,c) = customFilt(Window{i,c},Fs,[1 10],2);
        adchf(:,c) = abs(diff(chf(:,c))).^2;
    end
    WindowTags{i} = dataTags( start : start + winLen-1 );
    
    figure(2); hold on;
        title('Filtered EOG Data');
        plot(chf); hold off;
    figure(3);hold on; 
        plot(adchf); hold off;
    %%If surpasses threshold THEN classify and store data.
    thresholdCheck = find(max(chf)>THRESHOLD1);
    thresholdCheck2 = find(max(chf)>THRESHOLD2);
    getClass = [];
    if ~isempty(thresholdCheck) && isempty(thresholdCheck2)
        tagsPresent = unique(WindowTags{i},'stable')
        getClass = mode(WindowTags{i})
        if (getClass~=0 && getClass~=1 && getClass~=2)
            tYA(cnt,1) = getClass;
            F_ch1(cnt,:) = featureExtractionEOG2( chf(:,1)' );%%%% TODO: Make separate EOG Classifier/Feature Extraction. 
            F_ch2(cnt,:) = featureExtractionEOG2( chf(:,2)' );
            cnt = cnt + 1;
        end
    end
%     if ~isempty(thresholdCheck2)
%         tagsPresent = unique(WindowTags{i},'stable')
%         getClass = mode(WindowTags{i}); 
%         if(getClass==1 || getClass==2)
%             tYA(cnt,1) = getClass;
%             %%%% TODO: Make separate EOG Classifier/Feature Extraction. 
%             F_ch1(cnt,:) = featureExtractionEOG( chf(:,1)' );
%             F_ch2(cnt,:) = featureExtractionEOG( chf(:,2)' );
%             cnt = cnt + 1;
%         end
%     end
%{  
    if ~isempty(thresholdCheck) && isempty(thresholdCheck2)
        while isempty(getClass)
            commandwindow;
            gC = mode(WindowTags{i})
            getClass = input('Enter an integer value!\n');
            if isempty(getClass)
                getClass = mode(WindowTags{i});
                %if(getClass~=0 && getClass~=1 && getClass~=2)
                if (getClass~=0)
                    tagsPresent = unique(WindowTags{i},'stable')
                    tYA(cnt,1) = getClass;
                    F_ch1(cnt,:) = featureExtractionEOG2( chf(:,1)' );%%%% TODO: Make separate EOG Classifier/Feature Extraction. 
                    F_ch2(cnt,:) = featureExtractionEOG2( chf(:,2)' );
                    cnt = cnt + 1;
                end
                getClass = [];
            else
                break;
            end
        end
    end
%}

    clf(figNum);
    clf(3);
end

%horzcat: 36 features (12/ch)
tXA = [F_ch1 F_ch2];

clearvars -except tXA tYA
%% % % % TRAIN CLASSIFIER, KNN :: filename('knnclassification')
% % % % % Combine %{
tX = [];
tY = [];
if exist('marc_eyemove_tD2.mat','file')==2
    load('marc_eyemove_tD2');
end
tX = [tX;tXA];
tY = [tY;tYA];
AAA = [tX,tY];
save('marc_eyemove_tD2','tX','tY','AAA');
%}
%% PLOT:
numberfeatures=13;
for i=1:numberfeatures
    figure(i); hold on;
    plot(tX(:,i));
    plot(tX(:,i+13));
end

%% Load Training data for Classification Learner: 
%{
clear;clc;close all;
load('marc_eyemove_tD.mat');
tXtY = [tX tY];
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