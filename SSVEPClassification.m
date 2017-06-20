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
% test signals:
sti_f = [ 9.5, 15.15, 16.67, 18.52, 20.0 ];
t_length=4;                              % data length (4 s)
TW=1:1:t_length;
TW_p=round(TW*Fs);
N=2;    % number of harmonics
ref0=refsig(sti_f(1),Fs,t_length*Fs,N);
ref1=refsig(sti_f(2),Fs,t_length*Fs,N);
ref2=refsig(sti_f(3),Fs,t_length*Fs,N);
ref3=refsig(sti_f(4),Fs,t_length*Fs,N);
ref4=refsig(sti_f(5),Fs,t_length*Fs,N);
% %{
for i = 1:length(wStart)
    for j = 1:length(range)
        windowLength = range(j); 
        fprintf('From [%d] to [%d] \r\n',wStart(i),wStart(i)+(windowLength-1));
        [PredictedClass(j,i),PC2(j,i)] = classifySSVEP(X1(wStart(i):wStart(i)+(windowLength-1)), X2(wStart(i):wStart(i)+(windowLength-1)), PLOTDATA, THRESHOLD_FRACTION);
    end
    if(size(DATA,2) == 3)
        ActualClass(i) = mode(DATA(wStart(i):wStart(i)+(windowLength-1),3));
    end
    X12 = [X1(wStart(i):wStart(i)+(windowLength-1)), X2(wStart(i):wStart(i)+(windowLength-1))]';
    tw_length = windowLength/Fs;
    [wx0,wy0,r0]=cca(X12,ref0(:,1:TW_p(tw_length)));
    [wx1,wy1,r1]=cca(X12,ref1(:,1:TW_p(tw_length)));
    [wx2,wy2,r2]=cca(X12,ref2(:,1:TW_p(tw_length)));
    [wx3,wy3,r3]=cca(X12,ref3(:,1:TW_p(tw_length)));
    [wx4,wy4,r4]=cca(X12,ref4(:,1:TW_p(tw_length)));
    [v,idx]=max([max(r0),max(r1),max(r2),max(r3),max(r4)]);
    PC3(i) = idx-1;
end
figure(1); hold on;
if(size(DATA,2) == 3)
    plot(ActualClass);
    for j = 1:length(range)
        plot(PredictedClass(j,:));
        plot(PC2(j,:));
        plot(PC3(j,:));
        Compare = ActualClass == PredictedClass(j,:);% | ActualClass == PredictedClass2;
        Compare2 = ActualClass == PC2(j,:);
        Compare3 = ActualClass == PC3(j,:);
        Accuracy(j) = sum(Compare)/length(ActualClass);
        Accuracy2(j) = sum(Compare2)/length(ActualClass);
        Accuracy3(j) = sum(Compare3)/length(ActualClass);
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

