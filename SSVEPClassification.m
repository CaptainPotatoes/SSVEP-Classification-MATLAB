%% SSVEP CLASSIFICATION:
clear;clc;close all;
% LOAD TRAINING DATA: (tX, tY);
DATA = csvread('Matt_10Hz_null.csv');
Fs = 250;
% Generating Idealized Signals:
clear sig_ideal;
len = 5000; % 4s = testSignal(desiredF,ln,amplitude,Fs);
f = [0.0000, 10.0000,12.5000,15.1515,16.6667];
sig_ideal = zeros(4,len);
for i = 1:length(f)
    [sig_ideal(i,:),T] = testSignal(f(i),len);
    figure(5); hold on; plot(T,sig_ideal(i,:)); ylim([-1.2E-4 1.2E-4]);
    %TODO: Extract Features:
end
%% Feature Extraction: Expanding window method:
% Todo move to separate function 
close all;clc;
range = 250:60:2500;
F = [];
filtRange = [8 20];
for i = 1:5
    F = [F; featureExtractionSSVEP(sig_ideal(i,:),range,filtRange,0)];
end

%% Feature Extraction for Signal
X_samples = DATA(:,1);
F_s = featureExtractionSSVEP(X_samples, range, filtRange, [0]);

%% CCA
% load carbig;
% X = [Displacement Horsepower Weight Acceleration MPG];
% nans = sum(isnan(X),2) > 0;
% [A,B,r,U,V] = CCA(X(~nans,1:3),X(~nans,4:5));
% Compare to each
for i = 1:5
    start = 1+(i-1)*38;
    F0{i} = F(start:start+37,:);
    [~,~,~,U{i},V{i}] = CCA(F_s, F0{i});
    if~isempty(U{i})
        plot(U{i},V{i},'.r')
    end
end

% plot(U(:,1),V(:,1),'.')
