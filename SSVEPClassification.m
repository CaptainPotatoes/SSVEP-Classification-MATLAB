%% SSVEP CLASSIFICATION:
clear;clc;close all;
% [DATA,filename] = csvread('Subject1_SingleChannel_10Hz_to_16Hz.csv');
% [DATA,filename] = csvread('EEG_SSVEPData_2017.05.31_14.55.24.csv');
[DATA, filename] = csvread('Subject1_Trial1.1.csv');
Fs = 250;
X_1 = DATA(:,1);
X_2 = DATA(:,2);
h = 1/250;
t=0:h:(size(DATA,1)/250)-h;
range = 250:250:1000;
start = 1;
wStart = start:250:(length(X_1)-max(range));
PLOTDATA = 1==0;
% figure(6); hold on; 
% plot(t,DATA(:,3),'r'),ylabel('Class Label'),xlabel('Time (s)'),title('Target Class');
% Filt/Extract Spectrograms 
%{
for i = 1:length(wStart)
    start = (i-1)*115+1; 
    fin = start + 114;%S(start:fin,:)
    [T,F,STFT1{i}] = extractSpectrograms(X_1(wStart(i):wStart(i)+999),PLOTDATA);
    [~,~,STFT2{i}] = extractSpectrograms(X_2(wStart(i):wStart(i)+999),PLOTDATA);
    CLASS{i} = unique(DATA(wStart(i):wStart(i)+999,3));
%     a = input('Continue? \n');
end
clearvars -except filename T F STFT1 STFT2 CLASS
save([filename(1:end-4) '_spect.mat'],'-v7','T','F','STFT1','STFT2','CLASS')
%}
%% Feature Extraction for Signal
% filtRange = [8 20];
pts = [1, 7935, 15500, 23425];
start = pts(1);
% Generate table (reference):

wStart = start:250:(length(X_1)-max(range));
i=1;
% %%%%%%%%%%%%% % %{
PLOTDATA = 1==0;
THRESHOLD_FRACTION = 3;
for i = 1:length(wStart)
    CLASS(i) = classifySSVEP(X_1(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
end
for i = 1:length(wStart)
    CLASS2(i) = classifySSVEP(X_2(wStart(i):wStart(i)+999),PLOTDATA,THRESHOLD_FRACTION);
end
figure(7); 
h = 1/250;
t=0:h:(size(DATA,1)/250)-h;
hold on; 
plot(t,DATA(:,3),'r');
plot(CLASS),ylabel('Class Label'),xlabel('Time (s)'),title([filename '-Ch1']);
plot(CLASS2),legend('Target Class','Ch1','Ch2');
%}

%% CCA-KNN?
% Generate static reference signal
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
