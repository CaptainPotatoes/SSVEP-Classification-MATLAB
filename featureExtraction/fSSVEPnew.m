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
    %---FOR >= 500 DP ---%
threshPSDL = zeros(4,2);
threshPSDL(1,:) = [9.99 10.01];
threshPSDL(2,:) = [12.4 12.70]; 
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
fch1 = fch1(:);
fch2 = fch2(:);
fch3 = fch3(:);
windowLength = length(fch1);
fchw = zeros(nCh,windowLength);
fchw(1,1:windowLength) = fch1(1:windowLength);
fchw(2,1:windowLength) = fch2(1:windowLength);
fchw(3,1:windowLength) = fch3(1:windowLength);
% TODO: PREALLOCATE F BASED ON WINDOWLENGTH: (IF NECESSARY):

% all windows should be the same size:
FFT = zeros(4,1025);
select = false(4,1025);

wLFFT = zeros(4,1);
wLPSD = zeros(4,1);

% Data is already filtered:
if plotData
    fH = figure(1); %-% Figure Handle
    set(fH, 'Position', [2560, 0, 1280, 920]);
    xL = [9.0 17.2];
    clf(fH)
end

% TODO: FILL IN LATER!!!
%{%}
if windowLength < 500 && windowLength >= 250
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
                fft_sel_loc(ch,i) = fselect(fft_sel_L(1)); %verifies that maxes are peaks. [peak must occur w/i range]
                fft_sel_pks(ch,i) = fft_sel_P(1);
            else
                fft_sel_loc(ch,i) = 0;
                fft_sel_pks(ch,i) = 0;
            end
    %         
            if plotData
                subplot(2,2,1);hold on;
                plot(fselect, fftselect,'.b');
                plot(fftL(ch,i), fftM(ch,i), '^k');
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
        for i=1:4
            selectPSD(i,:) = fPSD>=threshPSD(i,1) & fPSD<=threshPSD(i,2);
            fPSDselect = fPSD(selectPSD(i,:));
            PSDselect = PSD(ch,selectPSD(i,:));
            [PSDM(ch,i), PSDL0(ch,i)] = max(PSDselect);
            PSDL(ch,i) = fPSDselect(PSDL0(ch,i));
            if length(PSDselect) > 3
                [psd_sel_P, psd_sel_L] = findpeaks(PSDselect,'SortStr','descend');
                if ~isempty(psd_sel_P)
                    psd_sel_loc(ch,i) = fPSD(psd_sel_L(1));
                    psd_sel_pks(ch,i) = psd_sel_P(1);
                else
                    psd_sel_loc(ch,i) = 0;
                    psd_sel_pks(ch,i) = 0;
                end
            else
                psd_sel_loc(ch,i) = 0;
                psd_sel_pks(ch,i) = 0;
            end
            if plotData
                subplot(2,2,2);hold on;
                plot(fPSDselect, PSDselect,'.b');
                plot(PSDL(ch,i), PSDM(ch,i), '^m');
                plot(psd_sel_loc(ch,i),psd_sel_pks(ch,i),'or');
            end
        end
        %TODO FIND PEAKS:
    end
    ch = 4;
    FFT(4,:) = (FFT(1,:)+FFT(2,:)+FFT(3,:));
    for i=1:4
            fselect = f(select(i,:));
            fftselect = FFT(ch,select(i,:));
            [fftM(ch,i), fftL0(ch,i)] = max(fftselect);
            fftL(ch,i) = fselect(fftL0(ch,i));
    %         
            [fft_sel_P, fft_sel_L] = findpeaks(fftselect,'SortStr','descend');
            if ~isempty(fft_sel_P)
                fft_sel_loc(ch,i) = fselect(fft_sel_L(1)); %verifies that maxes are peaks. [peak must occur w/i range]
                fft_sel_pks(ch,i) = fft_sel_P(1);
            else
                fft_sel_loc(ch,i) = 0;
                fft_sel_pks(ch,i) = 0;
            end
            if plotData
                subplot(2,2,1); hold on;
                plot(fft_sel_loc(ch,i),fft_sel_pks(ch,i),'or');
            end
    end
    if plotData
        subplot(2,2,1); hold on;
        plot(f, FFT(4,:)),xlim(xL);
    end
    PSD(4,:) = PSD(1,:)+PSD(2,:)+PSD(3,:);
    for i=1:4
        selectPSD(i,:) = fPSD>=threshPSD(i,1) & fPSD<=threshPSD(i,2);
        fPSDselect = fPSD(selectPSD(i,:));
        PSDselect = PSD(ch,selectPSD(i,:));
        [PSDM(ch,i), PSDL0(ch,i)] = max(PSDselect);
        PSDL(ch,i) = fPSDselect(PSDL0(ch,i));
        if length(PSDselect) > 3
            [psd_sel_P, psd_sel_L] = findpeaks(PSDselect,'SortStr','descend');
            if ~isempty(psd_sel_P)
                psd_sel_loc(ch,i) = fPSD(psd_sel_L(1));
                psd_sel_pks(ch,i) = psd_sel_P(1);
            else
                psd_sel_loc(ch,i) = 0;
                psd_sel_pks(ch,i) = 0;
            end
        else
            psd_sel_loc(ch,i) = 0;
            psd_sel_pks(ch,i) = 0;
        end
        if plotData
            subplot(2,2,2);hold on;
            plot(fPSD,PSD(ch,:)),xlim(xL);
            plot(fPSDselect, PSDselect,'.b');
            plot(PSDL(ch,i), PSDM(ch,i), '^m');
            plot(psd_sel_loc(ch,i),psd_sel_pks(ch,i),'or');
        end
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
                fft_sel_loc(ch,i) = fselect(fft_sel_L(1)); %verifies that maxes are peaks. [peak must occur w/i range]
                fft_sel_pks(ch,i) = fft_sel_P(1);
            else
                fft_sel_loc(ch,i) = 0;
                fft_sel_pks(ch,i) = 0;
            end
    %         
            if plotData
                subplot(2,2,1);hold on;
                plot(fselect, fftselect,'.b');
                plot(fftL(ch,i), fftM(ch,i), '^k');
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
        for i=1:4
            selectPSD(i,:) = fPSD>=threshPSDL(i,1) & fPSD<=threshPSDL(i,2);
            fPSDselect = fPSD(selectPSD(i,:));
            PSDselect = PSD(ch,selectPSD(i,:));
            if ~isempty(PSDselect)
                [PSDM(ch,i), PSDL0(ch,i)] = max(PSDselect);
                PSDL(ch,i) = fPSDselect(PSDL0(ch,i));
            else
                PSDM(ch,i) = 0;
                PSDL0(ch,i) = 0;
                PSDL(ch,i) = 0;
            end
            if length(PSDselect) > 3
                [psd_sel_P, psd_sel_L] = findpeaks(PSDselect,'SortStr','descend');
                if ~isempty(psd_sel_P)
                    psd_sel_loc(ch,i) = fPSD(psd_sel_L(1));
                    psd_sel_pks(ch,i) = psd_sel_P(1);
                else
                    psd_sel_loc(ch,i) = 0;
                    psd_sel_pks(ch,i) = 0;
                end
            else
                psd_sel_loc(ch,i) = 0;
                psd_sel_pks(ch,i) = 0;
            end
            if plotData
                subplot(2,2,2);hold on;
                plot(fPSDselect, PSDselect,'.b');
                plot(PSDL(ch,i), PSDM(ch,i), '^m');
                plot(psd_sel_loc(ch,i),psd_sel_pks(ch,i),'or');
            end
        end
        %TODO FIND PEAKS:
    end
    ch = 4;
    FFT(4,:) = (FFT(1,:)+FFT(2,:)+FFT(3,:));
    for i=1:4
            fselect = f(select(i,:));
            fftselect = FFT(ch,select(i,:));
            [fftM(ch,i), fftL0(ch,i)] = max(fftselect);
            fftL(ch,i) = fselect(fftL0(ch,i));
    %         
            [fft_sel_P, fft_sel_L] = findpeaks(fftselect,'SortStr','descend');
            if ~isempty(fft_sel_P)
                fft_sel_loc(ch,i) = fselect(fft_sel_L(1)); %verifies that maxes are peaks. [peak must occur w/i range]
                fft_sel_pks(ch,i) = fft_sel_P(1);
            else
                fft_sel_loc(ch,i) = 0;
                fft_sel_pks(ch,i) = 0;
            end
            if plotData
                subplot(2,2,1); hold on;
                plot(fft_sel_loc(ch,i),fft_sel_pks(ch,i),'or');
            end
    end
    if plotData
        subplot(2,2,1); hold on;
        plot(f, FFT(4,:)),xlim(xL);
    end
    PSD(4,:) = PSD(1,:)+PSD(2,:)+PSD(3,:);
    for i=1:4
%         selectPSD(i,:) = fPSD>=threshPSDL(i,1) & fPSD<=threshPSDL(i,2);
        fPSDselect = fPSD(selectPSD(i,:));
        PSDselect = PSD(ch,selectPSD(i,:));
        [PSDM(ch,i), PSDL0(ch,i)] = max(PSDselect);
        PSDL(ch,i) = fPSDselect(PSDL0(ch,i));
        if length(PSDselect) > 3
            [psd_sel_P, psd_sel_L] = findpeaks(PSDselect,'SortStr','descend');
            if ~isempty(psd_sel_P)
                psd_sel_loc(ch,i) = fPSD(psd_sel_L(1));
                psd_sel_pks(ch,i) = psd_sel_P(1);
            else
                psd_sel_loc(ch,i) = 0;
                psd_sel_pks(ch,i) = 0;
            end
        else
            psd_sel_loc(ch,i) = 0;
            psd_sel_pks(ch,i) = 0;
        end
        if plotData
            subplot(2,2,2);hold on;
            plot(fPSD,PSD(ch,:)),xlim(xL);
            plot(fPSDselect, PSDselect,'.b');
            plot(PSDL(ch,i), PSDM(ch,i), '^m');
            plot(psd_sel_loc(ch,i),psd_sel_pks(ch,i),'or');
        end
    end
    %Classification method #2 (w/ STFT):

    %TODO:
    h=64;
    wlen = 256;
    nfft = 2048;
    winLim = [9 17.6];
    K = sum(hammPeriodic(wlen))/wlen;
    [S1, F, T] = stft(fchw(1,:),wlen,h,nfft,Fs);
    [S2, ~, ~] = stft(fchw(2,:),wlen,h,nfft,Fs);
    [S3, ~, ~] = stft(fchw(3,:),wlen,h,nfft,Fs);
%     S1L = 20*log10(abs(S1)/wlen/K + 1E-6);
%     S2L = 20*log10(abs(S2)/wlen/K + 1E-6);
%     S3L = 20*log10(abs(S3)/wlen/K + 1E-6);
    SC = 20*log10(abs(S1(F<winLim(2) & F>winLim(1),:))/wlen/K + 1E-6)+ ...
        20*log10(abs(S2(F<winLim(2) & F>winLim(1),:))/wlen/K + 1E-6)+ ...
        20*log10(abs(S3(F<winLim(2) & F>winLim(1),:))/wlen/K + 1E-6);
    SummedRows = scaleAbs(sum(SC,2));
    F2 = F(F<winLim(2) & F>winLim(1));
    [M, I] = max(SummedRows);
    if plotData
        subplot(2,2,3)
        STFTlims = F<winLim(2) & F>winLim(1);
        imagesc(T,F2,SC),ylim(winLim);
        set(gca,'YDir','normal')
        handl = colorbar;
        colormap(jet);
        subplot(2,2,4);hold on;
        plot(F2, SummedRows);
        plot(F2(I), M, 'or');
    end
end
%     NEWFEATURES = ;
%% Collect Feature data into 'F'
    %First separate features by channel: (row vects)
    % first to remsove: *FFT_Ltop(2) ... not sure how I will use this
    % Also remove FFTPeaks2 and averageFFTPeak2
%WANT INFO TO PRINT IN ORDER:
    %% FPRINTFs:
if(plotData)
%     fprintf('Important Data (p2): [l = %d]\n',windowLength);

end
F = [1,2];
end %END FUNCTION

