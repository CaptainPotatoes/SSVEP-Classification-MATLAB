function [ P ] = fECONV2( X0, Fs, plotData )
X = zeros(1,length(X0));
X = (X0(:)');
ic = 0.1;
% Clusters First Set
C1 = 9.7:ic:10.3;
C2 = 12.1:ic:12.7;
C3 = 14.8:ic:15.4;
C4 = 16.2:ic:16.8;
% Second set (Updated).
% C1 = 8.8:ic:9.4;
% C2 = 10.8:ic:11.4;
% C3 = 16.4:ic:17;
% C4 = 18.2:ic:18.8;
% Setup (16.66; 18.5; 20; 5)
% C1 = 16.4:ic:17;
% C2 = 18.2:ic:18.8;
% C3 = 19.7:ic:20.3;
% C4 = 24.4:ic:25;

f_new = [C1,C2,C3,C4]; 
hannW = hannWin(2048);winLim = [6 30]; 
len = 2000;
fixed = zeros(1,len);
sigs = zeros(length(f_new),len);
convconv = zeros(length(f_new),length(X) + len - 1);
Mconv = zeros(1,length(f_new));
Lconv = Mconv;
if plotData
%     fH11 = figure(11); set(fH11, 'Position', [0, 0, 1440, 960]);
    fH = figure(12); %-% Figure Handle
    set(fH, 'Position', [0, 0, 1440, 960]);
    clf(fH)
end
% sigs = [];
for i = 1:length(f_new)
    [sigs(i,:)] = testSignal(f_new(i),len);
    fixed = sigs(i,:);
    convconv(i,:) = conv(X,fixed);
    [CPowerSpectrum ,wfreqs] = welch_psd(convconv(i,:), Fs, hannW); 
    [Mconv(1,i),Lconv(1,i)] = max( CPowerSpectrum ); 
    if plotData
        xlim(winLim); hold on; plot(wfreqs, CPowerSpectrum );
    end
end
if plotData
    hold on;plot(wfreqs(Lconv),Mconv,'*r'),xlim(winLim);
%     figure(11); plot(convconv');
end
P = [wfreqs(Lconv),Mconv];
end

