clear;clc;close all;
%      %%%%% FILE SELECTION %%%%%
DATA = csvread('Marc_Large_EyeMove_T1.csv');
% DATA = csvread('Marc_Trial_New_Blink1.csv');
        %%%%% - DATA SELECTION - %%%%%
%%%%%%-- THRESHOLDS --%%%%%%
    %%%_- FOR DIFF(FILT(CH))
UTH1 = 0.4E-4;
UTH2 = 3.25E-4;% UTH2 = 2.75E-4;
LTH1 = -0.5E-4;
LTH2 = -2.75E-4; %TODO: Change to lower number
rFB = 0; % Remove From Beginning
rFE = 0; % Remove From End
numCh = 2;
for i = 1:numCh
    ch(:,i) = DATA(1+rFB:end-rFE,i);
end
dataTags = DATA(1+rFB:end-rFE,3);
        %%%%% - VARIABLES & THRESHOLDS - %%%%%
TH1 = 2.85E-4; TH2 = 3.6E-3;TH3 = -0.5E-3;
DIFF_UPTH1 = 0.4E-4;
DIFF_LOTH1 = -0.5E-4;
Fs = 250; h=1/250; s = length(ch)/250; 
t=0:h:(s-h); % Time Signal
recordingLengthSeconds = length(ch(:,1))/250;
fprintf('Recording Length: %2.2f seconds \n',recordingLengthSeconds);
%%% PLOT RAW DATA RECORDING:
% figure(10); plot(ch),xlim([0,length(ch)]); 
%%% CALCULATE FILTERED AND DIFFERENCE ARRAYS %%%
for i = 1:numCh
%     filtch(:,i) = customFilt(ch(:,i),Fs,[0.15,9.5],3);
    filtch(:,i) = customFilt(ch(:,i),Fs,[1,9.9],5);
    diffchf(:,i) = diff(filtch(:,i));
end
%%%% PLOT FILTERED TIME SIGNAL %%%%
f1 = figure(1); set(f1, 'Position', [100, 100, 1600, 900]); hold on;
title('EOG Signal Amplitude');%
xlabel('Time (s)'); 
ylabel('EOG Signal Amplitude (mV)');

plot(t,filtch),xlim([0,max(t)]);hold on;legend('EOG Channel 1','EOG Channel 2')
for i=2:length(filtch)
    if (dataTags(i) ~= dataTags(i-1)) %Signal Change
        text( t(i), filtch(i,1), num2str(dataTags(i)) );
    end
end
%%%% PLOT DIFFERENTIAL TIME SIGNAL %%%
f2 = figure(2); set(f2, 'Position', [2600, 100, 1600, 900]); hold on; 
plot(diffchf);
rl1 = refline(0,LTH1); rl1.Color = 'r'; 
rl2 = refline(0,UTH1); rl2.Color = 'r';
rl3 = refline(0,LTH2); rl3.Color = 'm';
rl4 = refline(0,UTH2); rl4.Color = 'm';
for i=2:length(filtch)
    if (dataTags(i) ~= dataTags(i-1)) %Signal Change
        text( i, diffchf(i,1), num2str(dataTags(i)) );
    end
end
hold off;

%% Set classes:
    %%%%% - SET UP FIGURES - %%%%%
figNumA = 3;
figNumB = 4;
% fH = figure(figNumA); 
% set(fH, 'Position', [100, 100, 1600, 900]);
fH2 = figure(figNumB);
% set(fH2, 'Position', [2600, 100, 1600, 900]);
set(fH2, 'Position', [100, 100, 1600, 900]);
    %%%%% - Classification Vars - %%%%%
seconds = 4; %5 second window
winLen = seconds*Fs; 
winFraction = 4; %2.5; %1/4 of a second
winShift = floor( Fs/winFraction ); 
dataLimit = floor( (length(ch)-winLen)/winLen );
% Window = cell( seconds*winFraction*dataLimit - 1, numCh );
WindowTags = cell( 1,1 ); 
iterations = seconds*winFraction*dataLimit;
fprintf('Iterations: %2.2f \r\n',iterations);
chf = zeros(seconds*Fs,numCh);
selCh = chf;
dchf = zeros(seconds*Fs-1,numCh);
cnt = 1;

for i = 1 : iterations
    start = 1 + winShift*(i-1);
    winEnd = start + winLen-1;
    fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    for c = 1:numCh %         Window{i,c} = ch(start : start + winLen-1,c);
        selCh(:,c) = ch(start : start + winLen-1,c); %Selected Data
        chf(:,c) = customFilt(selCh(:,c),Fs,[0.15,9.5],3);
        dchf(:,c) = diff(chf(:,c)); % abs(diff(chf(:,c))).^2;
    end
    WindowTags{i} = dataTags(winEnd-249:winEnd);%WindowTags{i} = dataTags( start : start + winLen-1 );
%     figure(figNumA); hold on; title('Filtered EOG Data'); plot(chf(end-249:end,:)); hold off;
    figure(figNumB);hold on; title('Diff(Filtered EOG Data)'); plot(dchf(end-249:end,:)); 
    refline(0,UTH1); ar1 = refline(0,UTH2); refline(0,LTH1); ar2 = refline(0,LTH2); ar1.Color = 'r';ar2.Color='r';
    thresholdCheck = find(max(dchf(end-249:end,:))>UTH1);
    thresholdCheck2 = find(max(dchf(end-249:end,:))>UTH2);
    % Feature Extraction For Eye Movement
    F1(cnt,:) = featureExtractionEOG3( dchf(end-249:end,1), LTH1, LTH2, UTH1, UTH2, true );
    F2(cnt,:) = featureExtractionEOG3( dchf(end-249:end,2), LTH1, LTH2, UTH1, UTH2, true );
    % Feature Extraction For Eye Blinking
    F3(cnt,:) = featureExtractionEOG2( dchf(end-249:end,1), LTH1, LTH2, UTH1, UTH2, true );
    F4(cnt,:) = featureExtractionEOG2( dchf(end-249:end,2), LTH1, LTH2, UTH1, UTH2, true );
    getClass = [];
    % ALL CLASSES %{
% :: 4/26 :: TODO CHANGE SO EVERY WINDOW GETS A CLASS NUM!
    if isempty(thresholdCheck) && isempty(thresholdCheck2)
%         tYA(cnt,1) = 0;
%         cnt = cnt + 1;
    else
       tagsPresent = unique(WindowTags{i},'stable')
       getClass = mode(WindowTags{i})
       commandwindow; getInput = input('Enter an integer value!\n'); % Approve/Disapprove? 0 = reject:
       if isempty(getInput)
            tYA(cnt,1) = getClass;
            cnt = cnt + 1
       else
            if getInput~=-1
                tYA(cnt,1) = getInput;
                cnt = cnt + 1
            %else -1; rejected samples
            end
       end
    end
    %}
    % BLINK/DB CLASSES     
    %{
    if ~isempty(thresholdCheck) && ~isempty(thresholdCheck2)
        tagsPresent = unique(WindowTags{i},'stable')
        getClass = mode(WindowTags{i})
        commandwindow; getInput = input('Enter an integer value!\n'); % Approve/Disapprove?
        if isempty(getInput)
            if (getClass==1 || getClass==2)
                tYA(cnt,1) = getClass;
                cnt = cnt + 1;
            end
        end
    end
    %}
    % NULL CLASS
    %{
    if isempty(thresholdCheck) && isempty(thresholdCheck2)
        tagsPresent = unique(WindowTags{i},'stable')
        getClass = mode(WindowTags{i})
        if (getClass==0)
            tYA(cnt,1) = getClass;
            cnt = cnt + 1;
        end
    end
    %}
    % DIRECTION CLASSES     
    %{
    if ~isempty(thresholdCheck) && isempty(thresholdCheck2)
        tagsPresent = unique(WindowTags{i},'stable')
        getClass = mode(WindowTags{i})
        commandwindow; getInput = input('Enter an integer value!\n'); % Approve/Disapprove?
        if isempty(getInput)
            if (getClass~=0 && getClass~=1 && getClass~=2)
                tYA(cnt,1) = getClass;
                cnt = cnt + 1;
            end
        end
    end
    %}
%     clf(figNumA);
    clf(figNumB);
end
tX0 = [F1 F2 F3 F4]; %F1 F2 
tXA = tX0(1:length(tYA),:);
tXB = [F1 F2 F3 F4];
% tXC = [F1 F2 F3 F4 L1 L2];
clearvars -except tXA tYA tXB F1 F2 F3 F4
AAA=[tXA tYA];

%% % % % TRAIN CLASSIFIER, KNN :: filename('knnclassification')
% % % % % Combine 
%{
filename = 'marc_eyemove_tD427.mat';
tX = [];
tY = [];
if exist(filename ,'file')==2
    -
end
tX = [tX;tXA];
tY = [tY;tYA];
AAA = [tX,tY];
save(filename,'tX','tY','AAA');
%}
%% PLOT: 
%{

numberfeatures=13;
for i=1:numberfeatures
    figure(i); hold on;
    plot(tX(:,i));
    plot(tX(:,i+13));
end
%}
%% Load Training data for Classification Learner: 
%{
clear;clc;close all;
load('marc_eyemove_tD.mat');
tXtY = [tX tY];

AAB = AAA;
for i=1:size(AAB,1)
    if AAB(i,27) == 0
        AAB(i,:) = [];
    end
end
%Split AAA;
AAA = AAB;
AS = size(AAA,2);
tX = AAA(:,1:AS-1);
tY = AAA(:,AS);
%}