function [ F ] = fSSVEPnew( fch1, fch2, fch3, Fs, plotData )
% Feature Extraction Function for Tri-channel SSVEP Feature Extraction:
% ----- INPUTS -----
% fch1, fch2, fch3: Tri-channel SSVEP Samples of certain window size
% MUST BE A VECTOR!
% MUST BE FILTERED!!
% Fs - Sampling Rate. 
% Using 250-sample windows, feature extraction is obtained using FFT and
% PSD

    %%%%% - Thresholds: - %%%%%
%----FFT----%
threshFFT = zeros(4,2);
threshFFT(1,:) = [9.5 10.63];%-% windows around certain target frequencies
threshFFT(2,:) = [11.9 12.7]; 
threshFFT(3,:) = [14.6 15.5];
threshFFT(4,:) = [16.2 16.74];
wLFFT = zeros(4,1);
    %---FOR >= 500 DP ---%
threshFFTL = zeros(4,2);
threshFFTL(1,:) = [9.76  10.14];%-% windows around certain target frequencies
threshFFTL(2,:) = [12.2  12.57]; 
threshFFTL(3,:) = [14.76 15.27];
threshFFTL(4,:) = [16.45 16.80];
    %Also use wLFFT
%----PSD----%
threshPSD = zeros(4,2);
threshPSD(1,:) = [9.5 10.5];
threshPSD(2,:) = [12 13]; 
threshPSD(3,:) = [14 15.1];
threshPSD(4,:) = [16 17];
wLPSD = zeros(4,1);
    %---FOR >= 500 DP ---%
threshPSDL = zeros(4,2);
threshPSDL(1,:) = [9.99 10.01];
threshPSDL(2,:) = [12.4 12.6]; 
threshPSDL(3,:) = [14.9 15.22];
threshPSDL(4,:) = [16 17];
%---STFT---%
    %---FOR >= 500 DP ---%
threshSTFT = zeros(4,2);
threshSTFT(1,:) = [9.60 10.26];
threshSTFT(2,:) = [12.2 12.71]; 
threshSTFT(3,:) = [14.75 15.39];
threshSTFT(4,:) = [16.11 16.98];
    %---FOR >= 1250 DP ---%
threshSTFT5 = zeros(4,2);
threshSTFT5(1,:) = [9.76  10.03];
threshSTFT5(2,:) = [12.3  12.5]; 
threshSTFT5(3,:) = [14.85 15.25];
threshSTFT5(4,:) = [16.33 16.7];
%----PREALLOCATE----%
nCh = 3;
FFTPeaks1 = zeros(1,nCh+1);
FFTPeaks2 = zeros(1,nCh+1);
PSDPeaks1 = zeros(1,nCh+1);
fch1 = fch1(:);
fch2 = fch2(:);
fch3 = fch3(:);
if plotData
%     figure(2);clf(2);
%     hold on;
%     plot(fch1);
%     plot(fch2);
%     plot(fch3);
%     hold off;
end
windowLength = length(fch1);
fchw = zeros(nCh,windowLength);
fchw(1,1:windowLength) = fch1(1:windowLength);
fchw(2,1:windowLength) = fch2(1:windowLength);
fchw(3,1:windowLength) = fch3(1:windowLength);
% TODO: PREALLOCATE F BASED ON WINDOWLENGTH: (IF NECESSARY):

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
select = false(4,1025);
% Data is already filtered:
if plotData
    fH = figure(1); %-% Figure Handle
    set(fH, 'Position', [2560, 0, 1280, 920]);
    xL = [9.0 17.2];
    clf(fH)
end
%between 250?500dp

% TODO: FILL IN LATER!!!
%{%}
if windowLength < 500
    for ch=1:nCh
    [f, FFT(ch,:)] = get_nfft_data(fchw(ch,:), Fs, 2048);
        if plotData
            subplot(2,2,1);hold on;
            plot(f,FFT(ch,:)),xlim(xL);
        end
        for i=1:4
            select(i,:) = f>threshFFT(i,1) & f<threshFFT(i,2);
            fselect = f(select(i,:));
            fftselect = FFT(ch,select(i,:));
            [fftM(ch,i), fftL0(ch,i)] = max(fftselect);
            fftL(ch,i) = fselect(fftL0(ch,i));
    %         
            [fft_sel_P, fft_sel_L] = findpeaks(fftselect,'SortStr','descend');
            if ~isempty(fft_sel_P)
                fft_sel_loc(ch,i) = fselect(fft_sel_L); %verifies that maxes are peaks. [peak must occur w/i range]
                fft_sel_pks(ch,i) = fft_sel_P(1);
            else
                fft_sel_loc(ch,i) = 0;
                fft_sel_pks(ch,i) = 0;
            end
    %         
            if plotData
                subplot(2,2,1);hold on;
                plot(fselect, fftselect,'.b');
                plot(fftL(ch,i), fftM(ch,i), '^c');
                plot(fft_sel_loc(ch,i),fft_sel_pks(ch,i),'or');
            end
        end
        len = windowLength;
        hW = hannWin(len);
        [PSD(ch,:), fPSD] = welch_psd(fchw(ch,:), Fs, hW);%fin-start
        %TODO FIND PEAKS:
    end
    ch = 4;
    FFT(4,:) = (FFT(1,:)+FFT(2,:)+FFT(3,:));
    for i=1:4
    %         select(i,:) = f>threshFFT(i,1) & f<threshFFT(i,2);
            fselect = f(select(i,:));
            fftselect = FFT(ch,select(i,:));
            [fftM(ch,i), fftL0(ch,i)] = max(fftselect);
            fftL(ch,i) = fselect(fftL0(ch,i));
    %         
            [fft_sel_P, fft_sel_L] = findpeaks(fftselect,'SortStr','descend');
            if ~isempty(fft_sel_P)
                fft_sel_loc(ch,i) = fselect(fft_sel_L); %verifies that maxes are peaks. [peak must occur w/i range]
                fft_sel_pks(ch,i) = fft_sel_P(1);
            else
                fft_sel_loc(ch,i) = 0;
                fft_sel_pks(ch,i) = 0;
            end
            if plotData
                subplot(2,2,3); hold on;
                plot(fft_sel_loc(ch,i),fft_sel_pks(ch,i),'or');
            end
    end
    if plotData
        subplot(2,2,3); hold on;
        plot(f, FFT(4,:)),xlim(xL);
    end
else %--- Data >=500 dp ---% 
    for ch=1:nCh
    [f, FFT(ch,:)] = get_nfft_data(fchw(ch,:), Fs, 2048);
        if plotData
            subplot(2,2,1);hold on;
            plot(f,FFT(ch,:)),xlim(xL);
        end
        for i=1:4
            select(i,:) = f>threshFFTL(i,1) & f<threshFFTL(i,2);
            fselect = f(select(i,:));
            fftselect = FFT(ch,select(i,:));
            [fftM(ch,i), fftL0(ch,i)] = max(fftselect);
            fftL(ch,i) = fselect(fftL0(ch,i));
    %         
            [fft_sel_P, fft_sel_L] = findpeaks(fftselect,'SortStr','descend');
            if ~isempty(fft_sel_P)
                fft_sel_loc(ch,i) = fselect(fft_sel_L); %verifies that maxes are peaks. [peak must occur w/i range]
                fft_sel_pks(ch,i) = fft_sel_P(1);
            else
                fft_sel_loc(ch,i) = 0;
                fft_sel_pks(ch,i) = 0;
            end
    %         
            if plotData
                subplot(2,2,1);hold on;
                plot(fselect, fftselect,'.b');
                plot(fftL(ch,i), fftM(ch,i), '^c');
                plot(fft_sel_loc(ch,i),fft_sel_pks(ch,i),'or');
            end
        end
        len = windowLength;
        hW = hannWin(len);
        [PSD(ch,:), fPSD] = welch_psd(fchw(ch,:), Fs, hW);%fin-start
        if plotData
        	subplot(2,2,2);hold on;
            plot(fPSD,PSD(ch,:)),xlim(xL);
        end
        %TODO FIND PEAKS FOR PSD:
    end
    ch = 4;
    FFT(4,:) = (FFT(1,:)+FFT(2,:)+FFT(3,:));
    for i=1:4
    %         select(i,:) = f>threshFFT(i,1) & f<threshFFT(i,2);
            fselect = f(select(i,:));
            fftselect = FFT(ch,select(i,:));
            [fftM(ch,i), fftL0(ch,i)] = max(fftselect);
            fftL(ch,i) = fselect(fftL0(ch,i));
    %         
            [fft_sel_P, fft_sel_L] = findpeaks(fftselect,'SortStr','descend');
            if ~isempty(fft_sel_P)
                fft_sel_loc(ch,i) = fselect(fft_sel_L); %verifies that maxes are peaks. [peak must occur w/i range]
                fft_sel_pks(ch,i) = fft_sel_P(1);
            else
                fft_sel_loc(ch,i) = 0;
                fft_sel_pks(ch,i) = 0;
            end
            if plotData
                subplot(2,2,3); hold on;
                plot(fft_sel_loc(ch,i),fft_sel_pks(ch,i),'or');
            end
    end
    if plotData
        subplot(2,2,3); hold on;
        plot(f, FFT(4,:)),xlim(xL);
    end
end
%{
for ch=1:nCh
% #1 Take FFT:
    % #1.1 Find Peaks and M/I
    %TODO NEW METHOD HERE:
 
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
ch = 4;
%Combine data into 'fourth' channel:
FFT(4,:) = (FFT(1,:)+FFT(2,:)+FFT(3,:));
for i=1:4
%         select(i,:) = f>threshFFT(i,1) & f<threshFFT(i,2);
        fselect = f(select(i,:));
        fftselect = FFT(ch,select(i,:));
        [fftM(ch,i), fftL0(ch,i)] = max(fftselect);
        fftL(ch,i) = fselect(fftL0(ch,i));
%         
        [fft_sel_P, fft_sel_L] = findpeaks(fftselect,'SortStr','descend');
        if ~isempty(fft_sel_P)
            fft_sel_loc(ch,i) = fselect(fft_sel_L); %verifies that maxes are peaks. [peak must occur w/i range]
            fft_sel_pks(ch,i) = fft_sel_P(1);
        else
            fft_sel_loc(ch,i) = 0;
            fft_sel_pks(ch,i) = 0;
        end
        if plotData
            subplot(2,2,3); hold on;
            plot(fft_sel_loc(ch,i),fft_sel_pks(ch,i),'or');
        end
end

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

if plotData
    subplot(2,2,3); hold on;
    plot(f, FFT(4,:)),xlim(xL);
%     plot(f(FFT_L), FFT_PKS, 'o');
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

b3 = (wLPSD(1)~=0) && (wLPSD(2)~=0) && (wLPSD(3)~=0) && (wLPSD(4)~=0);
if b3 %if a signal was detected for PSD on all channels:
    %check they are equivalent:
    testVect2 = [ wLPSD(1) wLPSD(1) wLPSD(1) wLPSD(1) ];
    b4 = isequal([wLPSD(1), wLPSD(2), wLPSD(3), wLPSD(4)], testVect2);
else
    b4 = false;
end

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
%}
%     NEWFEATURES = ;
%% Collect Feature data into 'F'
    %First separate features by channel: (row vects)
    % first to remsove: *FFT_Ltop(2) ... not sure how I will use this
    % Also remove FFTPeaks2 and averageFFTPeak2
%WANT INFO TO PRINT IN ORDER:
    averagePkRatioFFT = mean(FFT_PkRatio);
    averagePkRatioPSD = mean(PSD_PkRatio);
    %% FPRINTFs:
if(plotData)
    fprintf('Important Data: [l = %d]\n',windowLength);
    fprintf('FFT Matching Class %d %d %d %d \n',wLFFT(1),wLFFT(2)...
        ,wLFFT(3),wLFFT(4));
    fprintf('PSD Matching Class %d %d %d %d \n',wLPSD(1),wLPSD(2)...
        ,wLPSD(3),wLPSD(4));
    fprintf('FFT Peak Ratio: %1.3f\n',FFT_PkRatio);
    fprintf('PSD Peak Ratio: %1.3f\n',PSD_PkRatio);
    fprintf('Avg FFTL: %1.3f \n',averageFFTPeak);
    fprintf('Avg PSDL: %1.3f \n',averagePSDPeak);
    fprintf('Booleans: [%d %d %d %d] \n',b1,b2,b3,b4);
    fprintf('Avg FFTPkRatio: %1.3f \n',averagePkRatioFFT);
    fprintf('Avg PSDPkRatio: %1.3f \n',averagePkRatioPSD);
end
% F = [wLFFT' wLPSD' FFT_PkRatio' PSD_PkRatio' averageFFTPeak averagePSDPeak FFTPeaks1 PSDPeaks1 b1 b2 b3 b4 ];
F = [1,2];
end %END FUNCTION

