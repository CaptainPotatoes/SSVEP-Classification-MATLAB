clear;clc;close all;
[DATA, filename] = csvread('Subject1_Trial3.1.csv');
% [DATA, filename] = csvread('BioRadio_Matt_15s.csv');
% [DATA, filename] = csvread('EEG_SSVEPData_2017.);
X1 = DATA(:,1);
X2 = DATA(:,2);
Fs = 250; h=1/Fs;
t = 0:h:(size(DATA,1)/250)-h;
% range = 500:250:1000;
range = 1000;
maxlen = max(range);
winHop = 125;
wStart = 1:winHop:(length(X1)-max(range));
PLOTDATA = 1==0
THRESHOLD_FRACTION = 2 ; i = 1;
% ic = 0.1; i=1;
% C1 = 14.8:ic:15.4;
% C2 = 16.4:ic:17;
% C3 = 18.2:ic:18.8;
% C4 = 19.7:ic:20.3;
% f_new = [C1,C2,C3,C4];
% sigs = generateTestSignal(f_new, 2000);
% %{
for i = 30:length(wStart)
    for j = 1:length(range)
        windowLength = range(j); 
        fprintf('From [%d] to [%d] \r\n',wStart(i),wStart(i)+(windowLength-1));
%         [PredictedClass(j,i),PC2(j,i)] = classifySSVEP2(X1(wStart(i):wStart(i)+(windowLength-1)), X2(wStart(i):wStart(i)+(windowLength-1)), PLOTDATA, THRESHOLD_FRACTION);
        [PredictedClass(j,i),PC2(j,i)] = classifySSVEP(X1(wStart(i):wStart(i)+(windowLength-1)), X2(wStart(i):wStart(i)+(windowLength-1)), PLOTDATA, THRESHOLD_FRACTION);
    end
    if(size(DATA,2) == 3)
        ActualClass(i) = mode(DATA(wStart(i):wStart(i)+(windowLength-1),3));
    end
end
figure(1); hold on;
if(size(DATA,2) == 3)
    plot(ActualClass);
    for j = 1:length(range)
        plot(PredictedClass(j,:));
%         plot(PC2(j,:));
        Compare = ActualClass == PredictedClass(j,:);% | ActualClass == PredictedClass2;
        Compare2 = ActualClass == PC2(j,:);
        Accuracy(j) = sum(Compare)/length(ActualClass);
        Accuracy2(j) = sum(Compare2)/length(ActualClass);
    end
    
    legend('Actual Class','Predicted Class'); xlabel('Classifier Calls (1/s)');
    ylabel('Class Label');
    title(['SSVEP Classifier Output; Accuracy: ' num2str(100*Accuracy(1)) '%']);
else
    plot(PredictedClass); plot(PC2);
    legend('Predicted Class', 'PC2'); xlabel('Classifier Calls (1/s)');
    ylabel('Class Label');
    title(['SSVEP Classifier Output; Accuracy: ??%']);
end

%}

