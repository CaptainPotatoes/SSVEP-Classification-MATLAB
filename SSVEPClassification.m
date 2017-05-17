%% SSVEP CLASSIFICATION:
clear;clc;close all;
% LOAD TRAINING DATA: (tX, tY);
DATA = csvread('Matt_16Hz_null.csv');
Fs = 250;
% Generating Idealized Signals:
clear sig_ideal;
len = 5000; % 4s = testSignal(desiredF,ln,amplitude,Fs);
f = [0.0000, 10.0000,12.5000,15.1515,16.6667];
sig_ideal = zeros(4,len);
X_samples = DATA(:,1);
for i = 1:length(f)
    [sig_ideal(i,:),T,sig_h1(i,:),sig_h2(i,:),SQW(i,:)] = testSignal(f(i),len);
    figure(5); hold on; plot(T,sig_ideal(i,:)); ylim([-1.2E-4 1.2E-4]);
    %TODO: Extract Features:
end
% Some signal manipulation:
for i=1:length(f)
    conv0(i,:) = conv(customFilt(X_samples,Fs,[8 20],3),sig_ideal(i,:));
    conv1(i,:) = conv(X_samples,sig_h1(i,:));
    conv2(i,:) = conv(X_samples,sig_h2(i,:));
%     conv1(i,:) = conv(X_samples,SQW(i,:));
%     TODO: TAKE AVERAGE P2P ACROSS SIGNAL (ALSO TRY PEAK2RMS)
    % Fail. doesn't work. Try to analyze instead.
    c0p2p(i) = peak2peak(conv0(i,1:1000));    
    c1p2p(i) = peak2peak(conv1(i,1:1000));
    c2p2p(i) = peak2peak(conv2(i,1:1000));
    figure;hold on;plot(conv1(i,:));
end
% Analysis:
clf(4);fH = figure(4);hold on; set(fH, 'Position', [0, 0, 1600, 900]);%Spect
hannW = hannWin(2048);winLim = [6 24];
for i = 1:5
    [S1,wfreqs] = welch_psd(conv0(i,:), Fs, hannW);     plot(wfreqs, S1),xlim(winLim);
end


%% Feature Extraction: Expanding window method:
% Todo move to separate function 
close all;clc;
range = 380:60:2500;
F = []; F1 = [];
filtRange = [8 20];
for i = 1:5
    F0{i} = featureExtractionSSVEP(sig_ideal(i,:),range,filtRange,0);
%     F1{i} = featureExtractionSSVEP(sig_h1(i,:),range,filtRange,0);
%     F2{i} = featureExtractionSSVEP(sig_h2(i,:),range,filtRange,0);
    FF{i} = [F0{i}(:,5:8),F0{i}(:,13:16),F0{i}(:,21:24)];
    FF{i}(:,end+1) = i-1;
%     F = [F; featureExtractionSSVEP(sig_ideal(i,:),range,filtRange,0)];
end
FALL = [FF{1};FF{2};FF{3};FF{4};FF{5}];

%% Feature Extraction for Signal
F_s = featureExtractionSSVEP(X_samples, range, filtRange, [0]);
F_s2 = [F_s(:,5:8),F_s(:,13:16),F_s(:,21:24)];
for i = 1:size(F_s,1)
    Y(i) = knn(F_s2(i,:),FALL(:,1:end-1),FALL(:,end),1);
end
%% CCA
% load carbig;
% X = [Displacement Horsepower Weight Acceleration MPG];
% nans = sum(isnan(X),2) > 0;
% [A,B,r,U,V] = CCA(X(~nans,1:3),X(~nans,4:5));
% Compare to each
for i = 1:5
    start = 1+(i-1)*38;
    [~,~,~,U{i},V{i}] = CCA(F_s2, F1{i});
    if~isempty(U{i})
        figure
        plot(U{i}(:,1),V{i}(:,1),'.r')
    end
end

% plot(U(:,1),V(:,1),'.')
