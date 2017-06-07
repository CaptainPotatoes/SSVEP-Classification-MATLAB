clear;clc;close all;
[DATA, filename] = csvread('EEG_16_18_20_25_2.csv');
% [DATA, filename] = csvread('data\Subject1_SingleChannel_10Hz_to_16Hz.csv'); % load('SSVEP_featuresA1.mat'); tX = tDY(:,1:40); tY = tDY(:,41);
% load('data1.mat'); filename = 'def';
% DATA = DATA(1:250*30,:);
X1 = DATA(:,1);
X2 = DATA(:,2);
Fs = 250; h=1/Fs;
t = 0:h:(size(DATA,1)/250)-h;
range = 250:250:1000;
maxlen = max(range);
winHop = 250;
wStart = 1:winHop:(length(X1)-max(range));
PLOTDATA = 1==0
freqs = [10,12.5,15.15,16.666];
THRESHOLD_FRACTION = 3; i = 1;
% ic = 0.1;
% C1 = 16.4:ic:17;
% C2 = 18.2:ic:18.8;
% C3 = 19.7:ic:20.3;
% C4 = 24.4:ic:25;
% f_new = [C1,C2,C3,C4];
% sigs = generateTestSignal(f_new, 2000);
% %{
for i = 1:length(wStart)
%     combined(i,:) = convoluteSignals(X1(wStart(i):wStart(i)+999), X2(wStart(i):wStart(i)+999));
    [PredictedClass(i)] = classifySSVEP3(X1(wStart(i):wStart(i)+999), X2(wStart(i):wStart(i)+999), PLOTDATA, 3.0);
    ActualClass(i) = mode(DATA(wStart(i):wStart(i)+999,3));
end
figure(1); plot(PredictedClass);% figure(2); plot(M1); 
% figure(1); hold on;  plot(PredictedClass2); figure(3); plot(M2);
Compare = ActualClass == PredictedClass;% | ActualClass == PredictedClass2;
Accuracy = sum(Compare)/length(ActualClass);
% NewThreshold = mean(mean(M1)); hold on;
% h = refline([0,NewThreshold]); h.Color = 'r';
%{
%     [PredictedClass1(i),FTS1(i,:), M1(i,:)] = classifySSVEP(X1(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
%     [PredictedClass2(i),FTS2(i,:), M2(i,:)] = classifySSVEP(X2(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
%}