function [ FEATURES ] = fESSVEP( fch1, fch2, fch3, fch4, Fs, plotData )
%Feature Extraction for SSVEP:

nCh = 4;
fch1 = fch1(:);
fch2 = fch2(:);
fch3 = fch3(:);
fch4 = fch4(:);
wL = length(fch1);
fch = zeros(nCh,wL); 
fch(1,1:wL) = fch1(1:wL);
fch(2,1:wL) = fch2(1:wL);
fch(3,1:wL) = fch3(1:wL);
fch(4,1:wL) = fch4(1:wL);
if plotData
    fH = figure(12); %-% Figure Handle
    set(fH, 'Position', [2560, 0, 1280, 1440]);
    xL = [9.0 17.2];
    clf(fH)
end
winLim = [9,18];
if mod(wL,2) == 1;
    PSD = zeros(4,(wL-1)/2);
else
    PSD = zeros(4,wL/2);
end
nfft = 2^nextpow2(wL);
FFT = zeros(4,(nfft/2)+1); %1025 for 2048
if wL >= 250 
    for ch=1:nCh
        [f, FFT(ch,:)] = get_nfft_data(fch(ch,:), Fs, wL);
        if plotData
            subplot(3,2,1);hold on;plot(f,FFT(ch,:)),xlim(xL);
        end
        hW = hannWin(wL);
        [PSD(ch,:), fPSD] = welch_psd(fch(ch,:), Fs, hW);%fin-start
        if plotData
            subplot(3,2,2);hold on;plot(fPSD,PSD(ch,:)),xlim(xL);
        end
    end
end
h=64;
wlen = 256;
F = (0:(ceil((1+nfft)/2))-1)*Fs/nfft;
select = F<winLim(2) & F>winLim(1);
F1 = zeros(1,sum(select));
c = size(0:h:(wL-wlen),2);
S1 = zeros(sum(select), c);
S2 = zeros(sum(select), c);
S3 = zeros(sum(select), c);
S4 = zeros(sum(select), c);
K = sum(hammPeriodic(wlen))/wlen;
if wL>=500
    %STFT Variables
    F1=F(select);
    %spect:
    [S, ~, T] = stft2(fch(1,:),wlen,h,nfft,Fs);
    S1 = 20*log10( abs( S(select,:) ) /wlen/K + 1E-6 );
    [S, ~, ~] = stft2(fch(2,:),wlen,h,nfft,Fs);
    S2 = 20*log10( abs( S(select,:) ) /wlen/K + 1E-6 );
    [S, ~, ~] = stft2(fch(3,:),wlen,h,nfft,Fs);
    S3 = 20*log10( abs( S(select,:) ) /wlen/K + 1E-6 );
    [S, ~, ~] = stft2(fch(3,:),wlen,h,nfft,Fs);
    S4 = 20*log10( abs( S(select,:) ) /wlen/K + 1E-6 );
    if plotData
        subplot(3,2,3);hold on;imagesc(T,F1,S1),ylim(winLim),xlim([min(T),max(T)]);set(gca,'YDir','normal');colorbar;colormap(jet);
        subplot(3,2,4);hold on;imagesc(T,F1,S2),ylim(winLim),xlim([min(T),max(T)]);set(gca,'YDir','normal');colorbar;colormap(jet);
        subplot(3,2,5);hold on;imagesc(T,F1,S3),ylim(winLim),xlim([min(T),max(T)]);set(gca,'YDir','normal');colorbar;colormap(jet);
        subplot(3,2,6);hold on;imagesc(T,F1,S4),ylim(winLim),xlim([min(T),max(T)]);set(gca,'YDir','normal');colorbar;colormap(jet);
    end
end
%% TEMP
% Feats = zeros(8,1);
FEATURES = [S1(1,1),S1(2,1),S1(3,1),S1(4,1),F1(1,1),F1(1,2),F1(1,3),F1(1,4)];
end %/function fESSVEP

