%% Coder Executor
% Executes all paths for generating C/C++ Code:
clear;close all;clc
%{
load('mssvep_t2_16_2');
A = Trial{1}(1:end-250,1); 
B = Trial{2}(1:end-250,1);
C = Trial{3}(1:end-250,1);
D = Trial{4}(1:end-250,1);
Fs = SamplingRate;
start = 1;
i=1;
for wL = 250:6:5000
    wL
    W1 = A(start:wL);
    W2 = B(start:wL);
    W3 = C(start:wL);
    W4 = D(start:wL);
    TY{i} = eegcfilt(W1);
%     Y_EOG{i} = fullHybridClassifier(W1,W2,W3,W4,Fs,true);
%     Y{i} = fullHybridClassifier(W1,W2,W3,W4,Fs,false);
    i = i+1;
end
%}
%%
%{
load carbig;
X = [Displacement Horsepower Weight Acceleration MPG];
nans = sum(isnan(X),2) > 0;
% [A,B,r,U,V] = canoncorr(X(~nans,1:3),X(~nans,4:5));
[A,B,r,U,V] = CCA(X(~nans,1:3),X(~nans,4:5));
plot(U(:,1),V(:,1),'.')
xlabel('0.0025*Disp+0.020*HP-0.000025*Wgt')
ylabel('-0.17*Accel-0.092*MPG')
%}

%% TODO: function [ Y ] = eogcfilt_a( X ) for ANDROID
ba = csvread('MarcTest0.csv');
bach1 = ba(1:1000,1);
bach1f = eogcfilt_a(bach1);
