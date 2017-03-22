%% Coder Executor
% Executes all paths for generating C/C++ Code:
clear;close all;clc

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
    Y_EOG{i} = fullHybridClassifier(W1,W2,W3,W4,Fs,true);
    Y{i} = fullHybridClassifier(W1,W2,W3,W4,Fs,false);
    i = i+1;
end
