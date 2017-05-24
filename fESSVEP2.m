function [ P ] = fESSVEP2( X0, plotData )
%FESSVEP Feature Extraction for single (m x 1) SSVEP EEG Data Vector
%   X (m x 1) vectorize input:
% Fix X size:
X = zeros(1,length(X0));
X = (X0(:)');
Fs = 250;
%%Convolution Amplification:
ic = 0.1; % Increment, % f_new=9:0.1:17;
% Clusters:
C1 = 9.7:ic:10.3;
C2 = 12.1:ic:12.7;
C3 = 14.8:ic:15.4;
C4 = 16.2:ic:16.8;
f_new = [C1,C2,C3,C4]; 
hannW = hannWin(2048);winLim = [6 24]; 
len = 2000; %5000, 2500, great, but slow. 2000 is quite fast
fixed = zeros(1,len);
sigs = zeros(length(f_new),len);
convconv = zeros(length(f_new),length(X) + len - 1);
Mconv = zeros(1,length(f_new));
Lconv = Mconv;
for i = 1:length(f_new)
    [sigs(i,:)] = testSignal(f_new(i),len);
end
for i = 1:length(f_new);
    convconv(i,:) = conv(X,sigs(i,:),'full');
    [CPowerSpectrum,wfreqs] = welch_psd(convconv(i,:), Fs, hannW); 
    [Mconv(1,i),Lconv(1,i)] = max(CPowerSpectrum ); 
    if plotData
        subplot(3,2,[5 6]),xlim(winLim); hold on; plot(wfreqs, CPowerSpectrum );
    end
end
if plotData
    subplot(3,2,[5 6]); hold on;plot(wfreqs(Lconv),Mconv,'*r'),xlim(winLim);
end
% Fts:
P = [wfreqs(Lconv),Mconv];
end

