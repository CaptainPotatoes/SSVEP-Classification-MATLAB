%% SSVEP CLASSIFICATION:
clear;clc;close all;
% LOAD TRAINING DATA: (tX, tY);
% [DATA,filename] = csvread('Subject1_SingleChannel_10Hz_to_16Hz.csv');
% [DATA,filename] = csvread('Subject1_Trial1.1.csv');
[DATA,filename] = csvread('Subject1_Trial1.3.csv');
Fs = 250;
X_1 = DATA(:,1);
X_2 = DATA(:,2);
figure(6); 
h = 1/250;
t=0:h:(size(DATA,1)/250)-h;
hold on; 
plot(t,DATA(:,3),'r'),ylabel('Class Label'),xlabel('Time (s)'),title('Target Class');
%% Feature Extraction for Signal
range = 250:250:1000;
% filtRange = [8 20];
% pts = [1, 7935, 15500, 23425];
% start = pts(1);
% Generate table (reference):
start = 1;
wStart = start:250:(length(X_1)-max(range));
i=1;
PLOTDATA = 1==0;
THRESHOLD_FRACTION = 3;
% CLASS(i) = classifySSVEP(X_samples(1:1000),1);
for i = 1:length(wStart)
    CLASS(i) = classifySSVEP(X_1(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
    CLASS2(i) = classifySSVEP(X_2(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
end
figure(6); 
h = 1/250;
t=0:h:(size(DATA,1)/250)-h;
hold on; 
plot(t,DATA(:,3),'r');
plot(CLASS),ylabel('Class Label'),xlabel('Time (s)'),title([filename '-Ch1']);
plot(CLASS2),legend('Target Class','Ch1','Ch2');
%% CCA-KNN?
%{
for i = 1:5
    start = 1+(i-1)*38;
    [~,~,~,U{i},V{i}] = CCA(F_s2, F1{i});
    if~isempty(U{i})
        figure
        plot(U{i}(:,1),V{i}(:,1),'.r')
    end
end
%}
% plot(U(:,1),V(:,1),'.')
% F_s2 = [F_s(:,5:8),F_s(:,13:16),F_s(:,21:24)];
% for i = 1:size(F_s,1)
%     Y(i) = knn(F_s2(i,:),FALL(:,1:end-1),FALL(:,end),1);
% end
