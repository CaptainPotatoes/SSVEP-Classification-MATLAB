clear;clc;close all;
% [DATA, filename] = csvread('Subject1_15_20_withAlpha_1.csv');
% [DATA, filename] = csvread('BioRadio_Matt_15s_2.csv');
[DATA, filename] = csvread('EEG_SSVEPData_2017.06.12_13.03.54.csv');
X1 = DATA(:,1);
X2 = DATA(:,2);
Fs = 250; h=1/Fs;
t = 0:h:(size(DATA,1)/250)-h;
range = 250:250:1000;
maxlen = max(range);
winHop = 250;
windowLength = 1000;
wStart = 1:winHop:(length(X1)-max(range));
PLOTDATA = 1==1
THRESHOLD_FRACTION = 1.5; i = 1;
% ic = 0.1; i=1;
% C1 = 14.8:ic:15.4;
% C2 = 16.4:ic:17;
% C3 = 18.2:ic:18.8;
% C4 = 19.7:ic:20.3;
% f_new = [C1,C2,C3,C4];
% sigs = generateTestSignal(f_new, 2000);
% %{
for i = 1:length(wStart)
    fprintf('From [%d] to [%d] \r\n',wStart(i),wStart(i)+(windowLength-1));
%     [PredictedClass(i)] = classifySSVEP5(X1(wStart(i):wStart(i)+999), X2(wStart(i):wStart(i)+999), PLOTDATA, THRESHOLD_FRACTION);
    [PredictedClass(i),PC(i)] = classifySSVEP2(X1(wStart(i):wStart(i)+(windowLength-1)), X2(wStart(i):wStart(i)+(windowLength-1)), PLOTDATA, THRESHOLD_FRACTION);
    if(size(DATA,2) == 3)
        ActualClass(i) = mode(DATA(wStart(i):wStart(i)+(windowLength-1),3));
    end
end
if(size(DATA,2) == 3)
    figure(1); plot(ActualClass); 
    Compare = ActualClass == PredictedClass;% | ActualClass == PredictedClass2;
    Accuracy = sum(Compare)/length(ActualClass);
end
figure(1); hold on; plot(PredictedClass);
legend('Predicted Class', 'Actual Class'); xlabel('Classifier Calls (1/s)');
ylabel('Class Label');
title(['SSVEP Classifier Output; Accuracy: ' num2str(100*Accuracy) '%']);
%}

