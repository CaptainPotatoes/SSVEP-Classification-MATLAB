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
threshFFT = zeros(4,2);
threshFFT(1,:) = [9.6 10.4];%-% windows around certain target frequencies
threshFFT(2,:) = [11.9 12.7]; 
threshFFT(3,:) = [14.6 15.5];
threshFFT(4,:) = [16.2 16.7];
wLFFT = zeros(4,1);
%----PSD----%
threshPSD = zeros(4,2);
threshPSD(1,:) = [9.9 10.1];
threshPSD(2,:) = [12 13]; 
threshPSD(3,:) = [14.9 15.1];
threshPSD(4,:) = [16 17];
wLPSD = zeros(4,1);
%----PREALLOCATE----%
nCh = 3;
FFTPeaks1 = zeros(1,nCh+1);
FFTPeaks2 = zeros(1,nCh+1);
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
%between 250?500dp
len = windowLength;
PSD = zeros(4,len/2);
if windowLength>=250 && windowLength<500
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
        % #2 Take PSD Estimate: (Welch method)
        % Prepare hanning window:
        
        hW = hannWin(len);
        [PSD(ch,:), fPSD] = welch_psd(fchw(ch,1:len), Fs, hW);%fin-start
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
%         PSD_L(4,:) = PSD_L(4,:)(:);
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
    for chn = 1:nCh+1
        FFTPeaks1(1,chn) = FFT_Ltop(chn,1);
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
        b3 = isequal([wLFFT(1), wLFFT(2), wLFFT(3), wLFFT(4)], testVect);
    else
        b3 = false;
    end
    averagePSDPeak = mean([PSD_Ltop(1,1) PSD_Ltop(2,1) ...
        PSD_Ltop(3,1) PSD_Ltop(4,1)]);
    %some alternate method???, 
    b2 = (wLPSD(1)~=0) && (wLPSD(2)~=0) && (wLPSD(3)~=0) && (wLPSD(4)~=0);
    if b2 %if a signal was detected for PSD on all channels:
        %check they are equivalent:
        testVect2 = [ wLPSD(1) wLPSD(1) wLPSD(1) wLPSD(1) ];
        b4 = isequal([wLPSD(1), wLPSD(2), wLPSD(3), wLPSD(4)], testVect2);
    else
        b4 = false;
    end
    
elseif windowLength>=500
    %Classification method #2 (w/ STFT):
    % Use CCA with longer time periods.
    %TODO:
else
    error('Not enough data!\n');
end %/windowLength>=250 && windowLength<500
    %% Collect Feature data into 'F'
    %First separate features by channel: (row vects)
    % first to remove: *FFT_Ltop(2) ... not sure how I will use this
    % Also remove FFTPeaks2 and averageFFTPeak2
    F0 = zeros(nCh+1, 8);
    for ch = 1:nCh+1
%         F(ch,:) = [FFT_Ltop(ch,1) FFT_Ltop(ch,2) FFT_MMM(ch,:) FFT_PkRatio(ch,:) ...
%             wLFFT(ch,:) PSD_Ltop(ch,1) PSD_Ltop(ch,2) PSD_MMM(ch,:) ...
%             PSD_PkRatio(ch,:) wLPSD(ch,:)];
        F0(ch,:) = [FFT_Ltop(ch,1) FFT_MMM(ch,:) FFT_PkRatio(ch,:) ...
            wLFFT(ch,:) PSD_Ltop(ch,1) PSD_MMM(ch,:) ...
            PSD_PkRatio(ch,:) wLPSD(ch,:)];
    end
%     Extras = [FFTPeaks1 FFTPeaks2 averageFFTPeak averageFFTPeak2 averagePSDPeak b1 b2];
    Extras = [FFTPeaks1 averageFFTPeak averagePSDPeak b1 b2 b3 b4];
    F = [F0(1,:) F0(2,:) F0(3,:) F0(4,:) Extras];
end %END FUNCTION
