%% Coder Executor
% Executes all paths for generating C/C++ Code:
clear;close all;clc
%{
load('mssvep_t2_16_2');
A = Trial{1}(1:end-250,1); 
B = Trial{2}(1:end-250,1);
C = Trial{3}(1:end-250,1);
D = Trial{4}(1:end-250,1);
Fs = SamplingRate;
start = 1;
i=1;
for wL = 250:6:5000
    wL
    W1 = A(start:wL);
    W2 = B(start:wL);
    W3 = C(start:wL);
    W4 = D(start:wL);
    TY{i} = eegcfilt(W1);
%     Y_EOG{i} = fullHybridClassifier(W1,W2,W3,W4,Fs,true);
%     Y{i} = fullHybridClassifier(W1,W2,W3,W4,Fs,false);
    i = i+1;
end
%}
%%
%{
load carbig;
X = [Displacement Horsepower Weight Acceleration MPG];
nans = sum(isnan(X),2) > 0;
% [A,B,r,U,V] = canoncorr(X(~nans,1:3),X(~nans,4:5));
[A,B,r,U,V] = CCA(X(~nans,1:3),X(~nans,4:5));
plot(U(:,1),V(:,1),'.')
xlabel('0.0025*Disp+0.020*HP-0.000025*Wgt')
ylabel('-0.17*Accel-0.092*MPG')
%}

%% function [ Y ] = eogcfilt_b( X ) for ANDROID 
%{
clear;clc;close all;
ba = csvread('MarcTest0.csv');
bach1 = ba(1:1000,1);
bach1f = customFilt(bach1,250,[0.15,9.5],3);
bach2f = eogcfilt_b(bach1);
figure(100);hold on;
plot(bach1f);plot(bach2f);
%}
%%  TODO: TEST EOG CLASSIFIER WITH FOR ANDROID
clear;clc;close all;
%      %%%%% FILE SELECTION %%%%%
% DATA = csvread('MarcDiff2ch_T5.csv');
% DATA = csvread('Marc_Large_EyeMove_T4.csv');
% DATA = csvread('Marc_Baseline_B_DB.csv');
DATA = csvread('Marc_Trial_New_largeGap_T1.csv');
        %%%%% - DATA SELECTION - %%%%%
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
    filtch(:,i) = customFilt(ch(:,i),Fs,[0.15,9.5],3);
    diffchf(:,i) = diff(filtch(:,i));
end
%%%% PLOT FILTERED TIME SIGNAL %%%%
f1 = figure(1); set(f1, 'Position', [100, 100, 1600, 900]); hold on;
title('EOG Signal Amplitude');%
xlabel('Time (s)'); 
ylabel('EOG Signal Amplitude (mV)');

plot(t,filtch),xlim([0,max(t)]);hold on;legend('EOG Channel 1','EOG Channel 2')
% rl1 = refline(0,TH1); rl1.Color = 'r'; % rl2 = refline(0,TH2); rl2.Color = 'm'; % rl3 = refline(0,TH3); rl3.Color = 'c';
for i=2:length(filtch)
    if (dataTags(i) ~= dataTags(i-1)) %Signal Change
        text( t(i), filtch(i,1), num2str(dataTags(i)) );
    end
end
%%%% PLOT DIFFERENTIAL TIME SIGNAL %%%
f2 = figure(2); set(f2, 'Position', [2600, 100, 1600, 900]); hold on; 
plot(diffchf);
rl4 = refline(0,DIFF_UPTH1); rl4.Color = 'r';
rl5 = refline(0,DIFF_LOTH1); rl5.Color = 'r';
for i=2:length(filtch)
    if (dataTags(i) ~= dataTags(i-1)) %Signal Change
        text( i, diffchf(i,1), num2str(dataTags(i)) );
    end
end
hold off;

%Set classes:
    %%%%% - SET UP FIGURES - %%%%%
figNumA = 3;
figNumB = 4;
fH = figure(figNumA); 
set(fH, 'Position', [100, 100, 1600, 900]);
fH2 = figure(figNumB);
set(fH2, 'Position', [2600, 100, 1600, 900]);
    %%%%% - Classification Vars - %%%%%
seconds = 4; %5 second window
winLen = seconds*Fs; 
winFraction = 10; %2.5; %1/4 of a second
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
%%%%%%-- THRESHOLDS --%%%%%%
    %%%_- FOR DIFF(FILT(CH))
UTH1 = 0.4E-4;
UTH2 = 2.75E-4;
LTH1 = -0.5E-4;
LTH2 = -2.75E-4;
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
%     figure(figNumB);hold on; title('Diff(Filtered EOG Data)'); plot(dchf(end-249:end,:)); 
    refline(0,UTH1); refline(0,UTH2); refline(0,LTH1); refline(0,LTH2); 
%     refline(0,DIFF_UPTH1); refline(0,DIFF_LOTH1); hold off; 
    %%TODO; Check threshold in last quartile.
%     thresholdCheck = find(max(dchf(end-249:end,:))>UTH1);
%     thresholdCheck2 = find(max(dchf(end-249:end,:))>UTH2);
%     F1(cnt,:) = featureExtractionEOG3( dchf(end-249:end,1), LTH1, LTH2, UTH1, UTH2, true );
%     F2(cnt,:) = featureExtractionEOG3( dchf(end-249:end,2), LTH1, LTH2, UTH1, UTH2, true );
    getClass = [];
    Y(i) = EOGClassifier(selCh(:,1),selCh(:,2));
    clf(figNumA);
    clf(figNumB);
end