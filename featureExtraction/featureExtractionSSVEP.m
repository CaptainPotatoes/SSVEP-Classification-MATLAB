function [ F ] = featureExtractionSSVEP( fch1, fch2, fch3, Fs )
% Feature Extraction Function for Tri-channel SSVEP Feature Extraction:
% ----- INPUTS -----
% fch1, fch2, fch3: Tri-channel SSVEP Samples of certain window size
% MUST BE A VECTOR!
% MUST BE FILTERED!!
% Fs - Sampling Rate. 
% Using 250-sample windows, feature extraction is obtained using FFT and
% PSD
%----FFT----%
wLFFT = zeros(4,2);
wLFFT(1,:) = [9.6 10.4];%-% windows around certain target frequencies
wLFFT(2,:) = [11.9 12.7]; 
wLFFT(3,:) = [14.6 15.5];
wLFFT(4,:) = [16.2 16.7];
%----PSD----%
wLPSD = zeros(4,2);
wLPSD(1,:) = [9.9 10.1];
wLPSD(2,:) = [12 13]; 
wLPSD(3,:) = [14.9 15.1];
wLPSD(4,:) = [16 17];
%----PREALLOCATE----%
nCh = 3;
% Ch = cell(nCh+1, 1);
FFTPeaks1 = zeros(1,nCh+1);
FFTPeaks2 = zeros(1,nCh+1);
windowLength = length(fch1);
fchw = zeros(nCh,windowLength);
fchw(1,:) = fch1;
fchw(2,:) = fch2;
fchw(3,:) = fch3;
% all windows should be the same size:

% Data is already filtered:
%between 250?500dp
if windowLength>=250 && windowLength<500
    for ch=1:nCh
    % #1 Take FFT:
        [f, FFT(ch,:)] = get_nfft_data(fchw(ch,:), Fs, 2048);
        % #1.1 Find Peaks and M/I
        [FFT_PKS(ch,:), FFT_L(ch,:)] = findpeaks(FFT(ch,:),'SortStr','descend');
        if length(FFT_PKS(ch,:))>1
            %Peak max minus min
            FFT_Ltop(ch,:) = f(FFT_L(ch,1:2));
            for w = 1:4
                if FFT_Ltop(ch,1)>wLFFT(w,1) && FFT_Ltop(ch,1)<wLFFT(w,2) 
                    FFT_MMM(ch,:) = FFT_PKS(ch,1) - FFT_PKS(ch,2);
                    FFT_PkRatio(ch,:) = FFT_PKS(ch,1)/FFT_PKS(ch,2);
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
        % Prepare static hann window: (may need to change to hann(250))
        hW = hann(windowLength);
        [PSD(ch,:), fPSD] = welch_psd(fchw(ch,:), Fs, hW);%fin-start
        % #2.2 Find Peaks and Max
        [PSD_PKS(ch,:), PSD_L(ch,:)] = findpeaks(PSD(ch,:),'SortStr','descend');
        if length(PSD_PKS(ch,:))>1
            PSD_Ltop(ch,:) = fPSD(PSD_L(ch,1:2));
            for w = 1:4
                if PSD_Ltop(ch,1)>=wLPSD(w,1) && PSD_Ltop(ch,1)<=wLPSD(w,2)
                    PSD_MMM(ch,:) = PSD_PKS(ch,1) - PSD_PKS(ch,2);
                    PSD_PkRatio(ch,:) = PSD_PKS(ch,1) / PSD_PKS(ch,2);
                    wLPSD(ch,:) = w;
                    break;
                else
                    PSD_MMM(ch,:) = 0;
                    PSD_PkRatio(ch,:) = 0;
                    wLPSD(ch,:) = 0;
                end
            end
        end
    end
    %Combine data into 'fourth' channel:
    FFT(4,:) = (FFT(1,:)+FFT(2,:)+FFT(3,:));
    [FFT_PKS(4,:), FFT_L(4,:)] = findpeaks(Ch{4}.FFT,'SortStr','descend');
    if length(FFT_PKS(4,:))>1
        FFT_Ltop(4,:) = f(FFT_L(4,1:2));
        for w = 1:4
            if FFT_Ltop(4,1)>wLFFT(w,1) && FFT_Ltop(4,1)<wLFFT(w,2) 
                FFT_MMM(4,:) = FFT_PKS(4,1) - FFT_PKS(4,2);
                FFT_PkRatio(4,:) = FFT_PKS(4,1)/FFT_PKS(4,2);
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
    [PSD_PKS(4,:), PSD_L(4,:)] = findpeaks(Ch{4}.PSD,'SortStr','descend');
    if length(PSD_PKS(4,:))>1
%         PSD_L(4,:) = PSD_L(4,:)(:);
        PSD_Ltop(4,:) = fPSD(PSD_L(4,1:2));
        for w = 1:4
            if PSD_Ltop(4,1)>=wLPSD(w,1) && PSD_Ltop(4,1)<=wLPSD(w,2)
                PSD_MMM(4,:) = PSD_PKS(4,1) - PSD_PKS(4,2);
                PSD_PkRatio(4,:) = PSD_PKS(4,1) / PSD_PKS(4,2);
                wLPSD(4,:) = w;
                break;
            else
                PSD_MMM(4,:) = 0;
                PSD_PkRatio(4,:) = 0;
                wLPSD(4,:) = 0;
            end
        end
    end
    for chn = 1:nCh+1
        FFTPeaks1(1,chn) = FFT_Ltop(chn,1);
        FFTPeaks2(1,chn) = FFT_Ltop(chn,2);
    end
    averageFFTPeak = mean([FFT_Ltop(1,1) FFT_Ltop(2,1) ...
        FFT_Ltop(3,1) FFT_Ltop(4,1)]);
    averageFFTPeak2 = mean([FFT_Ltop(1,2) FFT_Ltop(2,2) ...
        FFT_Ltop(3,2) FFT_Ltop(4,2)]);
    b1 = (wLFFT(1,:)~=0) && (wLFFT(2,:)~=0) && ...
        (wLFFT(3,:)~=0) && (wLFFT(4,:)~=0);
    averagePSDPeak = mean([PSD_Ltop(1,1) PSD_Ltop(2,1) ...
        PSD_Ltop(3,1) PSD_Ltop(4,1)]);
    b2 = (wLPSD(1,:)~=0) && (wLPSD(2,:)~=0) && (wLPSD(3,:)~=0) && (wLPSD(4,:)~=0);
elseif windowLength>=500
    %Classification method #2 (w/ STFT):
    %TODO:
else
    error('Not enough data!\n');
end %/windowLength>=250 && windowLength<500
    %% Collect Feature data into 'F'
    %First separate features by channel: (row vects)
    % first to remove: *FFT_Ltop(2) ... not sure how I will use this
    % Also remove FFTPeaks2 and averageFFTPeak2
%     F_1 = [Ch{1}.FFT_Ltop(1) Ch{1}.FFT_Ltop(2) Ch{1}.FFT_MMM Ch{1}.FFT_PkRatio Ch{1}.wLFFT ...
%         Ch{1}.PSD_Ltop(1) Ch{1}.PSD_Ltop(2) Ch{1}.PSD_MMM Ch{1}.PSD_PkRatio Ch{1}.wLPSD]; %10 features
    F = zeros(nCh+1, 10);
    for ch = 1:nCh+1
        F(ch,:) = [FFT_Ltop(ch,1) FFT_Ltop(ch,2) FFT_MMM(ch,:) FFT_PkRatio(ch,:) ...
            wLFFT(ch,:) PSD_Ltop(ch,1) PSD_Ltop(ch,2) PSD_MMM(ch,:) ...
            PSD_PkRatio(ch,:) wLPSD(ch,:)];
    end
    Extras = [FFTPeaks1 FFTPeaks2 averageFFTPeak averageFFTPeak2 averagePSDPeak b1 b2];
    F = [F(1,:) F(2,:) F(3,:) F(4,:) Extras];
end %END FUNCTION

