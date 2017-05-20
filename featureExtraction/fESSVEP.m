function [ SSVEP_FEATURES ] = fESSVEP( X, Fs, plotData )
%FESSVEP Feature Extraction for single (m x 1) SSVEP EEG Data Vector
%   X (m x 1) vectorize input:
% X         = X(:);
% Fs        = Sampling Frequency;
% plotData  = Plot Data

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
%%% - Constants - %%%
selc = ['.m';'.b';'.m';'.k']; %select dot color; 
nCh = 1;
winLim = [9,18];
% - Variables - %
if plotData
    fH = figure(12); %-% Figure Handle
    set(fH, 'Position', [0, 0, 1440, 960]);
    xL = [9.0 17.2];
    clf(fH)
end
wL = length(X);
if mod(wL,2) == 1;
    PSD = zeros(1,(wL-1)/2);
else
    PSD = zeros(1,wL/2);
end

Lfft = zeros(1,4);
Pfft = zeros(1,4);
Lpsd = zeros(1,4);
Ppsd = zeros(1,4);
Lstft = zeros(1,4);
Pstft = zeros(1,4);

nfft = 2^nextpow2(wL);
FFT = zeros(4,(nfft/2)+1);
hW = hannWin(wL);
% Lfft = zeros(4,1);
if wL >= 250 %Recall Data is Already Filtered From Calling Method
    [f, FFT] = get_nfft_data(X, Fs, wL);
    [PSD, fPSD] = welch_psd(X, Fs, hW);
    if plotData
        subplot(3,2,1);hold on;plot(f,FFT),xlim(xL);
        subplot(3,2,2);hold on;plot(fPSD,PSD),xlim(xL);
    end
    for i=1:4
        [fselect, fftselect, Lfft(i), Pfft(i)] = get_fft_features(f,FFT,threshFFT(i,:));
        [fselect2, psdselect, Lpsd(i), Ppsd(i)] = get_psd_features(fPSD,PSD,threshPSD(i,:));
        if plotData
            subplot(3,2,1);hold on;plot(fselect,fftselect,selc(i,:)); plot(Lfft(i),Pfft(i),'or');
            subplot(3,2,2);hold on;plot(fselect2,psdselect, selc(i,:)); plot(Lpsd(i),Ppsd(i),'or');
        end
    end
end
h = 32;
wlen = 128;
if wL >= 500
    wlen = 256;  
end
F = (0:(ceil((1+nfft)/2))-1)*Fs/nfft;
select = F<winLim(2) & F>winLim(1);
F1 = zeros(1,sum(select));
c = size(0:h:(wL-wlen),2);
SS = zeros(sum(select),1);
S1 = zeros(sum(select),c);
K = sum(hammPeriodic(wlen))/wlen;
M = zeros(nCh,4);
I = zeros(nCh,4);

if wL>=250 
    F1 = F(select);
    % STFT:
    [S,~,T] = stft2(X,wlen,h,nfft,Fs);
    S1 = 20*log10( abs( S(select,:) ) /wlen/K + 1E-6 );
    SS(:,1) = sum(S1,2)/size(S1,2);
    for i=1:4
        [fselect, stftselect, Lstft(i), Pstft(i), M(i), I(i)] = get_stft_features(F1,SS,threshFFT(i,:));
        if plotData
            subplot(3,2,4);hold on;plot(fselect,stftselect,selc(i,:));
            if Lstft(i)~=0
                plot(Lstft(i), Pstft(i), 'or');
            end
                plot(fselect(I(i)),M(i),'*k');
        end
    end
    if plotData
        subplot(3,2,3);hold on;imagesc(T,F1,S1),ylim(winLim),xlim([min(T),max(T)]);set(gca,'YDir','normal');colorbar;colormap(jet);
        subplot(3,2,4);hold on;plot(F1,SS(:));
    end
end
%%Convolution Amplification:
f_new=9:0.1:17;
hannW = hannWin(2048);winLim = [6 24]; 
len = 5000;
for i = 1:length(f_new)
    [sigs(i,:)] = testSignal(f_new(i),len);
    convconv(i,:) = conv(X,sigs(i,:),'full');
    [S1 ,wfreqs] = welch_psd(convconv(i,:), Fs, hannW); 
    [Mconv(i),L(i)] = max(S1); 
    if plotData
        subplot(3,2,[5 6]); hold on; plot(wfreqs, S1);plot(wfreqs(L),Mconv,'*r'),xlim(winLim);
    end
end

% Fts:
SSVEP_FEATURES = [Lfft,Pfft,Lpsd,Ppsd,Lstft,Pstft];
% SSVEP_FEATURES_SHORT = [Lpsd,Pfft,Ppsd,Pstft]; 
end

