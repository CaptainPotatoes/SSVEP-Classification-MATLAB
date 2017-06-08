%% SSVEP CLASSIFICATION:
clear;clc;close all;
% [DATA,filename] = csvread('Subject1_SingleChannel_10Hz_to_16Hz.csv');
% [DATA,filename] = csvread('EEG_SSVEPData_2017.05.31_14.55.24.csv');
[DATA, filename] = csvread('data\Subject1_SingleChannel_10Hz_to_16Hz.csv');
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
plot(t,DATA(:,3),'r'),ylabel('Class Label'),xlabel('Time (s)'),title('Target Class');
%% Feature Extraction for Signal
% filtRange = [8 20];
pts = [1, 7935, 15500, 23425];
start = pts(1);
% Generate table (reference);
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