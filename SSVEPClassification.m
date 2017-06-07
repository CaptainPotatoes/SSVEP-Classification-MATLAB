clear;clc;close all;
% [DATA, filename] = csvread('EEG_SSVEPData_2017.06.06_15.49.14.csv');
[DATA, filename] = csvread('Subject2_16_25_1.csv');
% [DATA, filename] = csvread('data\Subject1_SingleChannel_10Hz_to_16Hz.csv'); % load('SSVEP_featuresA1.mat'); tX = tDY(:,1:40); tY = tDY(:,41);
% load('data1.mat'); filename = 'def';
% DATA = DATA(1:250*30,:);
X1 = DATA(:,1);
X2 = DATA(:,2);
Fs = 250; h=1/Fs;
t = 0:h:(size(DATA,1)/250)-h;
range = 250:250:1000;
maxlen = max(range);
winHop = 120;
wStart = 1:winHop:(length(X1)-max(range));
PLOTDATA = 1==0
freqs = [10,12.5,15.15,16.666];
%{
Y = zeros(length(freqs),maxlen); Y1 = Y;
for i = 1:length(freqs)
    [Y(i,:), T, Y1(i,:)] = testSignal(freqs(i), maxlen);
    [C(i),FTS(i,:)] = classifySSVEP(Y(i,:),PLOTDATA,4);
end 
Y = Y'; 
%}
% %{
THRESHOLD_FRACTION = 3; i = 1;
% %{
for i = 1:length(wStart)
%     [PredictedClass1(i),FTS1(i,:), M1(i,:)] = classifySSVEP(X1(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
%     [PredictedClass2(i),FTS2(i,:), M2(i,:)] = classifySSVEP(X2(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
    [PredictedClass(i)] = classifySSVEP3(X1(wStart(i):wStart(i)+999), X2(wStart(i):wStart(i)+999), PLOTDATA, 3.0);
    ActualClass(i) = mode(DATA(wStart(i):wStart(i)+999,3));
end
figure(1); plot(PredictedClass);% figure(2); plot(M1); 
% figure(1); hold on;  plot(PredictedClass2); figure(3); plot(M2);
Compare = ActualClass == PredictedClass;% | ActualClass == PredictedClass2;
Accuracy = sum(Compare)/length(ActualClass);
% NewThreshold = mean(mean(M1)); hold on;
% h = refline([0,NewThreshold]); h.Color = 'r';
%}