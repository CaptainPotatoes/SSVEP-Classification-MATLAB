%% SSVEP CLASSIFICATION:
clear;clc;close all;
% LOAD TRAINING DATA: (tX, tY);
[DATA,filename] = csvread('EEGTrainingData_2017.05.22_12.14.39.csv');
% DATA = csvread('Matt_1ch_10_to_16_3.csv');
% DATA = csvread('EEGTrainingData_A1.csv');
Fs = 250;
% Generating Idealized Signals:
clear sig_ideal;
len = 5000; % 4s = testSignal(desiredF,ln,amplitude,Fs);
sig_ideal = zeros(4,len);
X_samples = DATA(:,1);
%% Feature Extraction for Signal
range = 250:60:1000;
filtRange = [8 20];
pts = [1, 7935, 15500, 23425];
% start = pts(1);
start = 1;
% Generate table (reference):
wStart = start:250:(length(X_samples)-1000);
i=1;
% CLASS(i) = classifySSVEP(X_samples(1:1000),1);
for i = 1:length(wStart)
    CLASS(i) = classifySSVEP(X_samples(wStart(i):wStart(i)+999),false);
end
figure(6);plot(CLASS),ylabel('Class Label'),xlabel('Time (s)'),title(filename);
%     CLASS(i) = classifySSVEP(X_samples,wStart(i),Fs);
%     [ FS(i,:),CLASS(i) ] = classifySSVEP(X_samples,wStart(i),Fs);
%     [F(i,:),CLASS(i)] = featureExtractionSSVEP(X_samples,range,filtRange,[0],wStart(i),Fs);

% FALL=[F,CLASS'];
%%
%{
starti = 1500; l0 = 1000;
endi = starti + l0-1;
f_new=9:0.1:17;
fH = figure(4); hold on; set(fH, 'Position', [0, 0, 1600, 900]);%Spect
hannW = hannWin(2048);winLim = [6 24]; winLim2 = [18 35];
for i=1:length(f_new)
    [sigs(i,:)] = testSignal(f_new(i),len);
    X_filt_snip = customFilt(X_samples(starti:endi),Fs,[8 20],3);
    convconv(i,:) = conv(X_filt_snip,sigs(i,:),'full');
    [S1 ,wfreqs] = welch_psd(convconv(i,:), Fs, hannW); figure(4); hold on; plot(wfreqs, S1)
    [M,L] = max(S1); plot(wfreqs(L),M,'*r'),xlim(winLim);
%     [sigs_h1(i,:)] = testSignal(2*f_new(i),len);
%     X_filt_snip_h1 = customFilt(X_samples(starti:endi),Fs,[8 50],3);
%     convconv_h1(i,:) = conv(X_filt_snip_h1,sigs_h1(i,:),'full');
end
%}

% F_s2 = [F_s(:,5:8),F_s(:,13:16),F_s(:,21:24)];
% for i = 1:size(F_s,1)
%     Y(i) = knn(F_s2(i,:),FALL(:,1:end-1),FALL(:,end),1);
% end
%% Some signal manipulation: 
%{
f = [0.0000, 10.0000, 12.5000, 15.0500, 16.6667];%TODO: use for loop 
for i=1:length(f)
    X_filt_all = customFilt(X_samples,Fs,[8 20],3);
%     conv0(i,:) = conv(X_filt_all(starti:endi),sig_ideal(i,:),'full');%+(rand(1)-0.5)*0.000015
    X_filt_snip = customFilt(X_samples(starti:endi),Fs,[8 20],3);
    conv0(i,:) = conv(X_filt_snip,sig_ideal(i,:),'full');
%     conv1(i,:) = conv(X_samples,sig_h1(i,:));
%     conv2(i,:) = conv(X_samples,sig_h2(i,:));
%     TODO: TAKE AVERAGE P2P ACROSS SIGNAL (ALSO TRY PEAK2RMS)
    c0p2p(i) = peak2peak(conv0(i,1:1000));
end
% Analysis:
% clf(4);

wlen = 1024; h=64; nfft = 4096;
K = sum(hamming(wlen, 'periodic'))/wlen;
for i = 1:5
    [S1 ,wfreqs] = welch_psd(conv0(i,:), Fs, hannW); figure(4); hold on; plot(wfreqs, S1),xlim(winLim);
    [f, P1] = get_fft_data(conv0(i,:)',Fs); figure(2);hold on; plot(f,P1),xlim(winLim);
    figure; [S1, f1, t1] = stft( conv0(i,:), wlen, h, nfft, Fs ); S2 = 20*log10(abs(S1(f1<winLim(2) & f1>winLim(1),:))/wlen/K + 1e-6); 
    imagesc(t1,f1(f1<winLim(2) & f1>winLim(1)),S2);set(gca,'YDir','normal');xlabel('Time, s');ylabel('Frequency, Hz');colormap(jet);cb = colorbar;ylabel(cb, 'Power (db)')
end
%}

%% CCA
% load carbig;
% X = [Displacement Horsepower Weight Acceleration MPG];
% nans = sum(isnan(X),2) > 0;
% [A,B,r,U,V] = CCA(X(~nans,1:3),X(~nans,4:5));
% Compare to each 
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
