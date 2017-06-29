clear;clc;close all;
% [DATA, filename] = csvread('Subject1_Trial5.1.csv');
% [DATA, filename] = csvread('BioRadio_Matt_15s.csv');
[DATA, filename] = csvread('WheelchairControl_SSVEP_SecondDemo_Data_SubjectMP.csv');
X1 = DATA(:,1);
X2 = DATA(:,2);
Fs = 250; h=1/Fs;
t = 0:h:(size(DATA,1)/250)-h;
range = 250:250:1000;
% range = 1000;
maxlen = max(range);
winHop = 60;
wStart = 1:winHop:(length(X1)-max(range));
PLOTDATA = 1==11
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
PredictedClass = zeros(length(range),length(wStart));
PC2 = PredictedClass; PC3 = PC2; ActualClass = PC3; 
% %{
for i = 1:length(wStart)
    for j = 1:length(range)
        fprintf('From [%d] to [%d] \r\n',wStart(i),wStart(i)+(range(j)-1));
        [PredictedClass(j,i),PC2(j,i), Ppsd(i,:,j)] = classifySSVEP(X1(wStart(i):wStart(i)+(range(j)-1)), X2(wStart(i):wStart(i)+(range(j)-1)), PLOTDATA, THRESHOLD_FRACTION);
        if(size(DATA,2) == 3)
            ActualClass(j,i) = mode(DATA(wStart(i):wStart(i)+(range(j)-1),3));
        end
        X12 = [X1(wStart(i):wStart(i)+(range(j)-1)), X2(wStart(i):wStart(i)+(range(j)-1))]';
        tw_length = range(j)/Fs;
        [wx0,wy0,r0]=cca(X12,ref0(:,1:TW_p(tw_length)));
        [wx1,wy1,r1]=cca(X12,ref1(:,1:TW_p(tw_length)));
        [wx2,wy2,r2]=cca(X12,ref2(:,1:TW_p(tw_length)));
        [wx3,wy3,r3]=cca(X12,ref3(:,1:TW_p(tw_length)));
        [wx4,wy4,r4]=cca(X12,ref4(:,1:TW_p(tw_length)));
        [v,idx]=max([max(r0),max(r1),max(r2),max(r3),max(r4)]);
        PC3(j,i) = idx-1;
    end
end
Accuracy = zeros(1,length(range)); Accuracy2 = Accuracy; Accuracy3 = Accuracy;

if(size(DATA,2) == 3)
    for j = 1:length(range)
        figure(j); hold on;
        plot(ActualClass(j,:));
        plot(PredictedClass(j,:));
        plot(PC2(j,:));
        plot(PC3(j,:),'.');
        Compare = ActualClass(j,:) == PredictedClass(j,:);% | ActualClass == PredictedClass2;
        Compare2 = ActualClass(j,:) == PC2(j,:);
        Compare3 = ActualClass(j,:) == PC3(j,:);
        Accuracy(1,j) = sum(Compare)/length(ActualClass(j,:));
        Accuracy(2,j) = sum(Compare2)/length(ActualClass(j,:));
        Accuracy(3,j) = sum(Compare3)/length(ActualClass(j,:));
        legend('Actual Class','Predicted Class','PC2','CCA PC3'); xlabel('Classifier Calls (1/s)');
        ylabel('Class Label');
        title(['SSVEP Classifier Output (', num2str(range(j)) ,' samples) ; Accuracy: ' num2str(100*max(Accuracy(:,j))) '%']);
    end
else
    plot(PredictedClass); plot(PC2);
    legend('Predicted Class', 'PC2'); xlabel('Classifier Calls (1/s)');
    ylabel('Class Label');
    title(['SSVEP Classifier Output; Accuracy: ??%']);
end
for i = 1:length(range)
    figure; %hold on;
    plot(Ppsd(:,1:end,i),'-*');
end
% C = confusionmat(ActualClass,PC2);

%}

