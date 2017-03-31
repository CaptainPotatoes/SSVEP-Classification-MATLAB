function [ FEATURES ] = fESSVEP( fch1, fch2, fch3, fch4, Fs, plotData )
%Feature Extraction for SSVEP:
% - Constants - %
%%%%% - Thresholds: - %%%%%
threshFFT = zeros(4,2);
threshFFT(1,:) = [9.5 10.63];
threshFFT(2,:) = [11.9 12.7]; 
threshFFT(3,:) = [14.6 15.5];
threshFFT(4,:) = [16.1 17.2];
threshPSD = zeros(4,2);
threshPSD(1,:) = [9.0 11.0];
threshPSD(2,:) = [11 13]; 
threshPSD(3,:) = [14 15.5];
threshPSD(4,:) = [15.5 17.5];
% - Variables - %
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
selc = ['.m';'.b';'.m';'.k']; %select dot color; 
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
Lfft = zeros(nCh,4);
Pfft = zeros(nCh,4);
Lpsd = zeros(nCh,4);
Ppsd = zeros(nCh,4);
Lstft = zeros(nCh,4);
Pstft = zeros(nCh,4);
nfft = 2^nextpow2(wL);
FFT = zeros(4,(nfft/2)+1);
hW = hannWin(wL);
if wL >= 250 
    for ch=1:nCh
        [f, FFT(ch,:)] = get_nfft_data(fch(ch,:), Fs, wL);
        [PSD(ch,:), fPSD] = welch_psd(fch(ch,:), Fs, hW);
        if plotData
            subplot(4,2,1);hold on;plot(f,FFT(ch,:)),xlim(xL);
            subplot(4,2,2);hold on;plot(fPSD,PSD(ch,:)),xlim(xL);
        end
        for i=1:4
            [fselect, fftselect, Lfft(ch,i), Pfft(ch,i)] = get_fft_features(f,FFT(ch,:),threshFFT(i,:));
            [fselect2, psdselect, Lpsd(ch,i), Ppsd(ch,i)] = get_psd_features(fPSD,PSD(ch,:),threshPSD(i,:));
            if plotData
                subplot(4,2,1);hold on;plot(fselect,fftselect,selc(i,:)); plot(Lfft(ch,i),Pfft(ch,i),'or');
                subplot(4,2,2);hold on;plot(fselect2,psdselect, selc(i,:)); plot(Lpsd(ch,i),Ppsd(ch,i),'or');
            end
        end
    end
end
h=64;
if wL<500
    wlen = 128;
else
    wlen = 256;
end
F = (0:(ceil((1+nfft)/2))-1)*Fs/nfft;
select = F<winLim(2) & F>winLim(1);
F1 = zeros(1,sum(select));
c = size(0:h:(wL-wlen),2);
S1 = zeros(sum(select),c);
SS = zeros(sum(select),4);
S2 = zeros(sum(select),c);
S3 = zeros(sum(select),c);
S4 = zeros(sum(select),c);
K = sum(hammPeriodic(wlen))/wlen;
M = zeros(nCh,4);
I = zeros(nCh,4);
if wL>=250
    %STFT Variables
    F1=F(select);
    %spect:
    [S, ~, T] = stft2(fch(1,:),wlen,h,nfft,Fs);
    S1 = 20*log10( abs( S(select,:) ) /wlen/K + 1E-6 );
    [S, ~, ~] = stft2(fch(2,:),wlen,h,nfft,Fs);
    S2 = 20*log10( abs( S(select,:) ) /wlen/K + 1E-6 );
    [S, ~, ~] = stft2(fch(3,:),wlen,h,nfft,Fs);
    S3 = 20*log10( abs( S(select,:) ) /wlen/K + 1E-6 );
    [S, ~, ~] = stft2(fch(4,:),wlen,h,nfft,Fs);
    S4 = 20*log10( abs( S(select,:) ) /wlen/K + 1E-6 );
    SS(:,1) = sum(S1,2)/size(S1,2);
    SS(:,2) = sum(S2,2)/size(S2,2);
    SS(:,3) = sum(S3,2)/size(S3,2);
    SS(:,4) = sum(S4,2)/size(S4,2);
    for ch=1:4
        for i=1:4
            [fselect, stftselect, Lstft(ch,i), Pstft(ch,i), M(ch,i), I(ch,i)] = get_stft_features(F1,SS(:,ch),threshFFT(i,:));
            if plotData
                subplot(4,2,[7 8]);hold on;plot(fselect,stftselect,selc(i,:));
                if Lstft(ch,i)~=0
                    plot(Lstft(ch,i), Pstft(ch,i), 'or');
                end
                    plot(fselect(I(ch,i)),M(ch,i),'*k');
            end
        end
    end
    if plotData
        subplot(4,2,3);hold on;imagesc(T,F1,S1),ylim(winLim),xlim([min(T),max(T)]);set(gca,'YDir','normal');colorbar;colormap(jet);
        subplot(4,2,4);hold on;imagesc(T,F1,S2),ylim(winLim),xlim([min(T),max(T)]);set(gca,'YDir','normal');colorbar;colormap(jet);
        subplot(4,2,5);hold on;imagesc(T,F1,S3),ylim(winLim),xlim([min(T),max(T)]);set(gca,'YDir','normal');colorbar;colormap(jet);
        subplot(4,2,6);hold on;imagesc(T,F1,S4),ylim(winLim),xlim([min(T),max(T)]);set(gca,'YDir','normal');colorbar;colormap(jet);
        subplot(4,2,[7 8]);hold on;plot(F1,SS(:,1));plot(F1,SS(:,2));plot(F1,SS(:,3));plot(F1,SS(:,4));
    end
end
%% TEMP
% Feats = zeros(8,1);
FEATURES = [Lfft(:) Pfft(:) Lpsd(:) Ppsd(:) Lstft(:) Pstft(:)];
end %/function fESSVEP

