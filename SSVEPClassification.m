clear;clc;close all;
[DATA, filename] = csvread('Subject1_15_20_withAlpha_2.csv');
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
ic = 0.1; i=1;
C1 = 14.8:ic:15.4;
C2 = 16.4:ic:17;
C3 = 18.2:ic:18.8;
C4 = 19.7:ic:20.3;
f_new = [C1,C2,C3,C4];
sigs = generateTestSignal(f_new, 2000);
% [PredictedClass(i)] = classifySSVEP4(X1(wStart(i):wStart(i)+999), X2(wStart(i):wStart(i)+999), PLOTDATA, THRESHOLD_FRACTION);
% %{
for i = 1:length(wStart)
%     [PredictedClass(i)] = classifySSVEP4(X1(wStart(i):wStart(i)+999), X2(wStart(i):wStart(i)+999), PLOTDATA, THRESHOLD_FRACTION);
    [PredictedClass(i)] = classifySSVEP5(X1(wStart(i):wStart(i)+999), X2(wStart(i):wStart(i)+999), PLOTDATA, THRESHOLD_FRACTION);
    ActualClass(i) = mode(DATA(wStart(i):wStart(i)+999,3));
end
figure(1); plot(PredictedClass); hold on; plot(ActualClass);
% figure(1); hold on;  plot(PredictedClass2); figure(3); plot(M2);
Compare = ActualClass == PredictedClass;% | ActualClass == PredictedClass2;
Accuracy = sum(Compare)/length(ActualClass);
legend('Predicted Class', 'Actual Class'); xlabel('Classifier Calls (1/s)');
ylabel('Class Label');
title(['SSVEP Classifier Output; Accuracy: ' num2str(100*Accuracy) '%']);

%}

