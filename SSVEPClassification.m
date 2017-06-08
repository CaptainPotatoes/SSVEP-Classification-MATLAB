clear;clc;close all;
% [DATA, filename] = csvread('Subject1_15s_2.csv');
% load('SSVEP_featuresA1.mat'); tX = tDY(:,1:40); tY = tDY(:,41);
[DATA, filename] = csvread('EEG_SSVEPData_2017.06.08_16.09.29.csv');
% DATA = DATA(1:end-250*4,:);
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
THRESHOLD_FRACTION = 4; i = 1;
ic = 0.1;
% C1 = 9.7:ic:10.3;
% C2 = 12.1:ic:12.7;
% C3 = 14.8:ic:15.4;
% C4 = 16.2:ic:16.8;
C1 = 14.8:ic:15.4;
C2 = 16.4:ic:17;
C3 = 18.2:ic:18.8;
C4 = 19.7:ic:20.3;

f_new = [C1,C2,C3,C4];
sigs = generateTestSignal(f_new, 2000);
% %{
for i = 1:length(wStart)
%     combined(i,:) = convoluteSignals(X1(wStart(i):wStart(i)+999), X2(wStart(i):wStart(i)+999));
%     [PredictedClass(i),FTS1(i,:), M1(i,:)] = classifySSVEP(X1(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
    [PredictedClass(i)] = classifySSVEP4(X1(wStart(i):wStart(i)+999), X2(wStart(i):wStart(i)+999), PLOTDATA, THRESHOLD_FRACTION);
    ActualClass(i) = mode(DATA(wStart(i):wStart(i)+999,3));
end
figure(1); plot(PredictedClass); hold on; plot(ActualClass);
% figure(1); hold on;  plot(PredictedClass2); figure(3); plot(M2);
Compare = ActualClass == PredictedClass;% | ActualClass == PredictedClass2;
Accuracy = sum(Compare)/length(ActualClass);
legend('Predicted Class', 'Actual Class'); xlabel('Classifier Calls (1/s)');
ylabel('Class Label');
title(['SSVEP Classifier Output; Accuracy: ' num2str(100*Accuracy) '%']);
%{
%     [PredictedClass1(i),FTS1(i,:), M1(i,:)] = classifySSVEP(X1(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
%     [PredictedClass2(i),FTS2(i,:), M2(i,:)] = classifySSVEP(X2(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
%}

