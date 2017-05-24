%% SSVEP CLASSIFICATION:
clear;clc;close all;
% LOAD TRAINING DATA: (tX, tY);
% [DATA,filename] = csvread('EEGTrainingData_2017.05.22_12.10.06.csv');
% [DATA,filename] = csvread('Matt_1ch_10_to_16_3.csv');
[DATA,filename] = csvread('A16HzOnly.csv');
Fs = 250;
X_samples = DATA(:,1);
%% Feature Extraction for Signal
range = 250:250:1000;
% filtRange = [8 20];
% pts = [1, 7935, 15500, 23425];
% start = pts(1);
% Generate table (reference):
start = 1;
wStart = start:250:(length(X_samples)-max(range));
i=1;
PLOTDATA = 1==0;
THRESHOLD_FRACTION = 2;
% CLASS(i) = classifySSVEP(X_samples(1:1000),1);
for i = 1:length(wStart)
    CLASS(i) = classifySSVEP(X_samples(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
end
figure(6);plot(CLASS),ylabel('Class Label'),xlabel('Time (s)'),title(filename);
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
