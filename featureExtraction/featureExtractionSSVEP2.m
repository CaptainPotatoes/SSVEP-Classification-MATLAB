function [ F ] = featureExtractionSSVEP2( fch1, fch2, fch3, Fs, plotData )
% Feature Extraction Function for Tri-channel SSVEP Feature Extraction:
% ----- INPUTS -----
% fch1, fch2, fch3: Tri-channel SSVEP Samples of certain window size
% MUST BE A VECTOR!
% MUST BE FILTERED!!
% Fs - Sampling Rate. 
% Using 250-sample windows, feature extraction is obtained using FFT and
% PSD
%----FFT----%
threshFFT = zeros(4,2);
threshFFT(1,:) = [9.5 10.63];%-% windows around certain target frequencies
threshFFT(2,:) = [11.9 12.7]; 
threshFFT(3,:) = [14.6 15.5];
threshFFT(4,:) = [16.2 16.74];
wLFFT = zeros(4,1);
%----PSD----%
threshPSD = zeros(4,2);
threshPSD(1,:) = [9.5 10.5];
threshPSD(2,:) = [12 13]; 
threshPSD(3,:) = [14.9 15.1];
threshPSD(4,:) = [16 17];
wLPSD = zeros(4,1);
%----PREALLOCATE----%
nCh = 3;
FFTPeaks1 = zeros(1,nCh+1);
FFTPeaks2 = zeros(1,nCh+1);
PSDPeaks1 = zeros(1,nCh+1);
fch1 = fch1(:);
fch2 = fch2(:);
fch3 = fch3(:);
windowLength = length(fch1);
fchw = zeros(nCh,windowLength);
fchw(1,1:windowLength) = fch1(1:windowLength);
fchw(2,1:windowLength) = fch2(1:windowLength);
fchw(3,1:windowLength) = fch3(1:windowLength);
% all windows should be the same size:
FFT = zeros(4,1025);
FFT_Ltop = zeros(4,2);
FFT_MMM = zeros(4,1);
FFT_PkRatio = zeros(4,1);
PSD = zeros(4,windowLength/2);
PSD_Ltop = zeros(4,2);
PSD_MMM = zeros(4,1);
PSD_PkRatio = zeros(4,1);
averageFFTPeak = 0;
averageFFTPeak2 = 0;
averagePSDPeak = 0;
b1 = false; %0?
b2 = false;
b3 = false;
b4 = false;
% Data is already filtered:
if plotData
    fH = figure(1); %-% Figure Handle
    set(fH, 'Position', [2560, 0, 1280, 920]);
    xL = [9.0 17.2];
    fprintf('Important Data: [l = %d]\n',windowLength);
    clf(fH)
end
%between 250?500dp

% Preallocate for spd:
for ch=1:nCh
% #1 Take FFT:
    [f, FFT(ch,:)] = get_nfft_data(fchw(ch,:), Fs, 2048);
    % #1.1 Find Peaks and M/I
    [FFT_PKS, FFT_L] = findpeaks(FFT(ch,:),'SortStr','descend');
    if length(FFT_PKS)>1
        %Peak max minus min
        FFT_Ltop(ch,:) = f(FFT_L(1:2));
        for w = 1:4
            if FFT_Ltop(ch,1)>threshFFT(w,1) && FFT_Ltop(ch,1)<threshFFT(w,2) 
                FFT_MMM(ch,:) = FFT_PKS(1) - FFT_PKS(2);
                FFT_PkRatio(ch,:) = FFT_PKS(1)/FFT_PKS(2);
                wLFFT(ch,:) = w;
                break;
            else
                FFT_MMM(ch,:) = 0;
                FFT_PkRatio(ch,:) = 0;
                wLFFT(ch,:) = 0;
            end
        end
    end
    if plotData
        subplot(2,2,1);hold on;
        plot(f,FFT(ch,:)),xlim(xL);
        plot(f(FFT_L), FFT_PKS, 'or');
    end
    % #2 Take PSD Estimate: (Welch method)
    % Prepare hanning window:
    len = windowLength;
    hW = hannWin(len);
    [PSD(ch,:), fPSD] = welch_psd(fchw(ch,:), Fs, hW);%fin-start
    % #2.2 Find Peaks and Max
    [PSD_PKS, PSD_L] = findpeaks(PSD(ch,:),'SortStr','descend');
    if length(PSD_PKS)>1
        PSD_Ltop(ch,:) = fPSD(PSD_L(1:2));
        for w = 1:4
            if PSD_Ltop(ch,1)>=threshPSD(w,1) && PSD_Ltop(ch,1)<=threshPSD(w,2)
                PSD_MMM(ch,:) = PSD_PKS(1) - PSD_PKS(2);
                PSD_PkRatio(ch,:) = PSD_PKS(1) / PSD_PKS(2);
                wLPSD(ch,:) = w;
                break;
            else
                PSD_MMM(ch,:) = 0;
                PSD_PkRatio(ch,:) = 0;
                wLPSD(ch,:) = 0;
            end
        end
    end
    if plotData
        subplot(2,2,2);hold on;
        plot(fPSD,PSD(ch,:)),xlim(xL);
        plot(fPSD(PSD_L), PSD_PKS, '^r');
    end
end

%Combine data into 'fourth' channel:
FFT(4,:) = (FFT(1,:)+FFT(2,:)+FFT(3,:));
[FFT_PKS, FFT_L] = findpeaks(FFT(4,:),'SortStr','descend');
if length(FFT_PKS)>1
    FFT_Ltop(4,:) = f(FFT_L(1:2));
    for w = 1:4
        if FFT_Ltop(4,1)>threshFFT(w,1) && FFT_Ltop(4,1)<threshFFT(w,2) 
            FFT_MMM(4,:) = FFT_PKS(1) - FFT_PKS(2);
            FFT_PkRatio(4,:) = FFT_PKS(1)/FFT_PKS(2);
            wLFFT(4,:) = w;
            break;
        else
            FFT_MMM(4,:) = 0;
            FFT_PkRatio(4,:) = 0;
            wLFFT(4,:) = 0;
        end
    end
end
PSD(4,:) = PSD(1,:)+PSD(2,:)+PSD(3,:);
[PSD_PKS, PSD_L] = findpeaks(PSD(4,:),'SortStr','descend');
if length(PSD_PKS)>1
    PSD_Ltop(4,:) = fPSD(PSD_L(1:2));
    for w = 1:4
        if PSD_Ltop(4,1)>=threshPSD(w,1) && PSD_Ltop(4,1)<=threshPSD(w,2)
            PSD_MMM(4,:) = PSD_PKS(1) - PSD_PKS(2);
            PSD_PkRatio(4,:) = PSD_PKS(1) / PSD_PKS(2);
            wLPSD(4,:) = w;
            break;
        else
            PSD_MMM(4,:) = 0;
            PSD_PkRatio(4,:) = 0;
            wLPSD(4,:) = 0;
        end
    end
end
fprintf('FFT Matching Class %d %d %d %d \n',wLFFT(1),wLFFT(2)...
    ,wLFFT(3),wLFFT(4));
fprintf('PSD Matching Class %d %d %d %d \n',wLPSD(1),wLPSD(2)...
    ,wLPSD(3),wLPSD(4));
fprintf('FFT Peak Ratio: %1.3f\n',FFT_PkRatio);
fprintf('PSD Peak Ratio: %1.3f\n',PSD_PkRatio);
if plotData
    subplot(2,2,3); hold on;
    plot(f, FFT(4,:)),xlim(xL);
    plot(f(FFT_L), FFT_PKS, 'o');
%     subplot(2,2,4); hold on;
%     plot(fPSD, PSD(4,:)),xlim(xL);
%     plot(fPSD(PSD_L), PSD_PKS, '^');
end
for chn = 1:nCh+1
    FFTPeaks1(1,chn) = FFT_Ltop(chn,1);
    PSDPeaks1(1,chn) = PSD_Ltop(chn,1);
    FFTPeaks2(1,chn) = FFT_Ltop(chn,2);
end
averageFFTPeak = mean([FFT_Ltop(1,1) FFT_Ltop(2,1) ...
    FFT_Ltop(3,1) FFT_Ltop(4,1)]);
fprintf('Avg FFTL: %1.3f \n',averageFFTPeak);
averageFFTPeak2 = mean([FFT_Ltop(1,2) FFT_Ltop(2,2) ...
    FFT_Ltop(3,2) FFT_Ltop(4,2)]);
b1 = (wLFFT(1)~=0) && (wLFFT(2)~=0) && ...
    (wLFFT(3)~=0) && (wLFFT(4)~=0);
if b1 %if a signal was detected on each FFT
    %check that they are all the same;
    testVect = [ wLFFT(1) wLFFT(1) wLFFT(1) wLFFT(1) ];
    b2 = isequal([wLFFT(1), wLFFT(2), wLFFT(3), wLFFT(4)], testVect);
else
    b2 = false;
end
averagePSDPeak = mean([PSD_Ltop(1,1) PSD_Ltop(2,1) ...
    PSD_Ltop(3,1) PSD_Ltop(4,1)]);
fprintf('Avg PSDL: %1.3f \n',averagePSDPeak);
b3 = (wLPSD(1)~=0) && (wLPSD(2)~=0) && (wLPSD(3)~=0) && (wLPSD(4)~=0);
if b3 %if a signal was detected for PSD on all channels:
    %check they are equivalent:
    testVect2 = [ wLPSD(1) wLPSD(1) wLPSD(1) wLPSD(1) ];
    b4 = isequal([wLPSD(1), wLPSD(2), wLPSD(3), wLPSD(4)], testVect2);
else
    b4 = false;
end
fprintf('Booleans: [%d %d %d %d] \n',b1,b2,b3,b4);
if windowLength>=500
    %Classification method #2 (w/ STFT):
    % Use CCA with longer time periods.(?)
    %TODO:
    h=64;
    wlen = 256;
    nfft = 2048;
    K = sum(hammPeriodic(wlen))/wlen;
    [S1, F, T] = stft(fchw(1,:),wlen,h,nfft,Fs);
    [S2, ~, ~] = stft(fchw(2,:),wlen,h,nfft,Fs);
    [S3, ~, ~] = stft(fchw(3,:),wlen,h,nfft,Fs);
    S1L = 20*log10(abs(S1)/wlen/K + 1E-6);
    S2L = 20*log10(abs(S2)/wlen/K + 1E-6);
    S3L = 20*log10(abs(S3)/wlen/K + 1E-6);
    winLim = [9 17.6];
    SC = 20*log10(abs(S1(F<winLim(2) & F>winLim(1),:))/wlen/K + 1E-6)+ ...
        20*log10(abs(S2(F<winLim(2) & F>winLim(1),:))/wlen/K + 1E-6)+ ...
        20*log10(abs(S3(F<winLim(2) & F>winLim(1),:))/wlen/K + 1E-6);
    SummedRows = scaleAbs(sum(SC,2));
    F2 = F(F<winLim(2) & F>winLim(1));
    [M, I] = max(SummedRows);
    if plotData
    subplot(2,2,4);hold on;
    plot(F2, SummedRows);
    plot(F2(I), M, 'or');
    end
    %TODO: ADD TO FEATURES AND USE IN FINAL DECISION!
end %/windowLength>=500
    
%% Collect Feature data into 'F'
    %First separate features by channel: (row vects)
    % first to remove: *FFT_Ltop(2) ... not sure how I will use this
    % Also remove FFTPeaks2 and averageFFTPeak2
%WANT INFO TO PRINT IN ORDER:
    averagePkRatioFFT = mean(FFT_PkRatio);
    averagePkRatioPSD = mean(PSD_PkRatio);
    fprintf('Avg FFTPkRatio: %1.3f \n',averagePkRatioFFT);
    fprintf('Avg PSDPkRatio: %1.3f \n',averagePkRatioPSD);
    F = [wLFFT' wLPSD' FFT_PkRatio' PSD_PkRatio' averageFFTPeak averagePSDPeak FFTPeaks1 PSDPeaks1 b1 b2 b3 b4 ];

end %END FUNCTION

