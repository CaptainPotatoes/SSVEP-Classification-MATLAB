clear;clc;close all;
% [DATA, filename] = csvread('');
% DATA = (DATA(1:250*15*3,:))
[DATA, filename] = csvread('C:\Users\mmahmood31\OneDrive - Georgia Institute of Technology\Publishing\_SSVEP\_SSVEPDataShared\S3\Subject3_Trial1.0.csv');
X1 = DATA(:,1);
X2 = DATA(:,2);
Fs = 250; h=1/Fs;
t = 0:h:(size(DATA,1)/250)-h;
% range = 250:250:1000;
range = 250;
range_cca = 250;
maxlen = max(range);
winHop = 60;
i = 1;
% while (i<(length(X1)-range))
%     
% end
wStart = 1:winHop:(length(X1)-max(range));
PLOTDATA = 1==0
THRESHOLD_FRACTION = 2 ; i = 1;
% Test Signals:
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
% frequencies = (0:(998/2-1))*250/998;
rAll = [250, 500, 1000, 2000, 4000];
rAllLen = rAll.*2;
for i = 1:length(wStart)
    i
    for j = 1:length(range)
        fprintf('From [%d] to [%d] \r\n',wStart(i),wStart(i)+(range(j)-1));
%         [Y(i,:), Ppsd(i,:,j), featurePoint(i,:)] = classifySSVEP2(X1(wStart(i):wStart(i)+(range(j)-1)), X2(wStart(i):wStart(i)+(range(j)-1)), PLOTDATA, THRESHOLD_FRACTION);
        
%         for p = 1:length(rAll)
%             [PSD{i,p}] = extractPowerSpectrum(X1(wStart(i):wStart(i)+(rAllLen(p)-1)), X2(wStart(i):wStart(i)+(rAllLen(p)-1)), rAll(p));
%         end
        [Y(i,:), PSDfull(i,:)] = classifySSVEP(X1(wStart(i):wStart(i)+(range(j)-1)), X2(wStart(i):wStart(i)+(range(j)-1)), PLOTDATA, THRESHOLD_FRACTION, 1);
%         [P, F] = extractPowerSpectrum2ch(X1(wStart(i):wStart(i)+(range(j)-1)), X2(wStart(i):wStart(i)+(range(j)-1)));
%         if(PLOTDATA)
%             figure(1); 
%             plot(f,P(i,:));
%         end
        if(size(DATA,2) == 3)
            ActualClass(j,i) = mode(DATA(wStart(i):wStart(i)+(range(j)-1),3));
        end
        X12 = [X1(wStart(i):wStart(i)+(range_cca(j)-1)), X2(wStart(i):wStart(i)+(range_cca(j)-1))]';
        tw_length = range_cca(j)/Fs;
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
        plot(Y(:,1),'.');
        plot(Y(:,2),'*');
        plot(PC3(j,:),'.');
        Compare1 = ActualClass(j,:) == Y(:,1)';% | ActualClass == PredictedClass2;
        Compare2 = ActualClass(j,:) == Y(:,2)';
        Compare3 = ActualClass(j,:) == PC3(j,:);
        Accuracy(1,j) = sum(Compare1)/length(ActualClass(j,:));
        Accuracy(2,j) = sum(Compare2)/length(ActualClass(j,:));
        Accuracy(3,j) = sum(Compare3)/length(ActualClass(j,:));
        legend('Actual Class','Y1','Y2','CCA PC3'); xlabel('Classifier Calls (1/s)');
        ylabel('Class Label');
%         title(['SSVEP Classifier Output (', num2str(range(j)) ,' samples) ; Accuracy: ' num2str(100*max(Accuracy(:,j))) '%']);
        title(['SSVEP Classifier Output (', num2str(range(j)) ,' samples) ; Accuracy: ' num2str(100*Accuracy(1)) '%']);
    end
else
    plot(PredictedClass); plot(PC2);
    legend('Predicted Class', 'PC2'); xlabel('Classifier Calls (1/s)');
    ylabel('Class Label');
    title(['SSVEP Classifier Output; Accuracy: ??%']);
end
Accuracy
%% PLot all data:
%{
figure(2); %clf(2);
correct = 'o';
incorrect = 'x';
classColors = {'r','g','b','m','k'};
hold on;
for i = 1: length(wStart)
    if(Compare1(i))
        plot(featurePoint(i,1), featurePoint(i,2), [correct,classColors{Y(i,1)+1}]);
    else
        plot(featurePoint(i,1), featurePoint(i,2), [incorrect,classColors{ActualClass(i)+1}]);
    end
end
title('Feature Distribution of Sample Data');
ylabel('Spectrum Max Power (W/Hz)');
xlabel('Frequency (Hz)');
% fP2 = featurePoint(featurePoint(:,1)>14.4,2);
% MFP2 = max(fP2);
% ylim([0 MFP2]);
%%%
%{
figure(3); %clf(3);
correct = 'o';
incorrect = 'x';
classColors = {'r','g','b','m','k'};
hold on;
for i = 1: length(wStart)
    if(Compare1(i))
        plot(featurePoint(i,1), featurePoint(i,2), [correct,classColors{Y(i,2)+1}]);
    else
        plot(featurePoint(i,1), featurePoint(i,2), [incorrect,classColors{ActualClass(i)+1}]);
    end
end
title('Feature Distribution of Sample Data');
ylabel('Spectrum Max Power (W/Hz)');
xlabel('Frequency (Hz)');
xlim([14.0 21]);
%}

% xlim[
%
C = confusionmat(ActualClass,Y(:,1));
C3(1,:,:) = C;
ATN = sum(squeeze(sum(C3,1)),2);%sum(sum(C));
SumDim1 = squeeze(sum(C3,1));
figure;
for i=1:5
    PercentC(i,:) = (SumDim1(i,:)./(ATN(i)))*100;
end
labels = {'Alpha','15','16','18','20'};
heatmap(PercentC, labels, labels, '%0.2f%%','Colormap','jet','ShowAllTicks',0,'UseLogColorMap',true,'Colorbar',true,'ColorLevels',30,'MaxColorValue',100,'MinColorValue',0);
figure; 
heatmap(PercentC, labels, labels, [],'Colormap','jet','ShowAllTicks',0,'UseLogColorMap',true,'Colorbar',true,'ColorLevels',30,'MaxColorValue',100,'MinColorValue',0);
% OverallAccuracy = sum(diag(PercentC))/5
xlabel('Predicted'),ylabel('Actual');
title(['Confusion Matrix for conv-PSDA method, Accuracy = ', num2str(Accuracy(1)*100), '%']); 
% Accuracy
%}

