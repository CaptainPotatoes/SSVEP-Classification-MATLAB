function [ Ppsd ] = fESSVEP( X0, Fs, plotData )
%FESSVEP Feature Extraction for single (m x 1) SSVEP EEG Data Vector
%   X (m x 1) vectorize input:
% Fix X size:
X = zeros(1,length(X0));
X = (X0(:)');
% Fs        = Sampling Frequency;
% plotData  = Plot Data

%%%%% - Thresholds: - %%%%%
NUMBER_CLASSES = 5;
threshFFT = zeros(NUMBER_CLASSES,2);
threshFFT(1,:) = [9.0 12];
threshFFT(2,:) = [14.6 15.7]; 
threshFFT(3,:) = [16.1 17.2];
threshFFT(4,:) = [18.0 19.2];
threshFFT(5,:) = [19.3 20.7];
threshPSD = zeros(NUMBER_CLASSES,2);
threshPSD(1,:) = [9 13];
threshPSD(2,:) = [14 15.5];
threshPSD(3,:) = [16.1 17.2];
threshPSD(4,:) = [18.0 18.8];
threshPSD(5,:) = [19.4 20.2];
%%% - Constants - %%%
selc = ['.m';'.b';'.m';'.k';'.c']; %select dot color; 
nCh = 1;
winLim = [6,24];
% - Variables - %
if plotData
    fH = figure(12); %-% Figure Handle
    set(fH, 'Position', [0, 0, 1440, 960]);
    clf(fH)
end
wL = length(X);
if mod(wL,2) == 1;
    PSD = zeros(1,(wL-1)/2);
else
    PSD = zeros(1,wL/2);
end

Lfft = zeros(1,NUMBER_CLASSES);
Pfft = zeros(1,NUMBER_CLASSES);
Lpsd = zeros(1,NUMBER_CLASSES);
Ppsd = zeros(1,NUMBER_CLASSES);
Lstft = zeros(1,NUMBER_CLASSES);
Pstft = zeros(1,NUMBER_CLASSES);

nfft = 2^nextpow2(wL);
FFT = zeros(1,(nfft/2)+1);
hW = hannWin(wL);
if wL >= 250
    [f, FFT] = get_nfft_data(X, Fs, wL);
    [PSD, fPSD] = welch_psd(X, Fs, hW);
    if plotData
        subplot(3,2,1);hold on;plot(f,FFT),xlim(winLim);
        subplot(3,2,2);hold on;plot(fPSD,PSD),xlim(winLim);
    end
    for i=1:NUMBER_CLASSES
        [fselect, fftselect, Lfft(i), Pfft(i)] = get_fft_features(f,FFT,threshFFT(i,:));
        [fselect2, psdselect, Lpsd(i), Ppsd(i)] = get_psd_features(fPSD,PSD,threshPSD(i,:));
        if plotData
            subplot(3,2,1);hold on;plot(fselect,fftselect,selc(i,:)); plot(Lfft(i),Pfft(i),'or'); title('FFT Analysis');
            subplot(3,2,2);hold on;plot(fselect2,psdselect, selc(i,:)); plot(Lpsd(i),Ppsd(i),'or'); title('Power Spectral Density Est.');
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
M = zeros(nCh,NUMBER_CLASSES);
I = zeros(nCh,NUMBER_CLASSES);

if wL>498 
    F1 = F(select);
    % STFT:
    [S,~,T] = stft(X,wlen,h,nfft,Fs);
    S1 = 20*log10( abs( S(select,:) ) /wlen/K + 1E-6 );
    SS(:,1) = sum(S1,2)/size(S1,2);
    for i=1:NUMBER_CLASSES
        [fselect, stftselect, M(i), I(i)] = get_stft_features(F1,SS,threshFFT(i,:));
        if plotData
            subplot(3,2,4);hold on;plot(fselect,stftselect,selc(i,:));
%             if Lstft(i)~=0
%                 plot(Lstft(i), Pst`ft(i), 'or');
%             end
            if(I(i)~=0)
                plot(fselect(I(i)),M(i),'or');
            end
        end
    end
    if plotData
        subplot(3,2,3);hold on;imagesc(T,F1,S1),ylim(winLim),xlim([min(T),max(T)]);set(gca,'YDir','normal');colorbar;colormap(jet); title('Spectrogram');
        subplot(3,2,4);hold on;plot(F1,SS(:)); title('Spectrogram Average Power');
        figure(13); hold on; plot(Lpsd, Ppsd, '*');
    end
end

%%Convolution Amplification:
% P = zeros(1,70);
% SSVEP_FEATURES = [Pfft,Ppsd,Pstft];

end

