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
        [f, Ch.FFT(ch,:)] = get_nfft_data(fchw(ch,:), Fs, 2048);
        % #1.1 Find Peaks and M/I
        [Ch.FFT_PKS(ch,:), Ch.FFT_L(ch,:)] = findpeaks(Ch.FFT(ch,:),'SortStr','descend');
        if length(Ch.FFT_PKS(ch,:))>1
            %Peak max minus min
            Ch.FFT_Ltop(ch,:) = f(Ch.FFT_L(ch,1:2));
            for w = 1:4
                if Ch.FFT_Ltop(ch,1)>wLFFT(w,1) && Ch.FFT_Ltop(ch,1)<wLFFT(w,2) 
                    Ch.FFT_MMM(ch,:) = Ch.FFT_PKS(ch,1) - Ch.FFT_PKS(ch,2);
                    Ch.FFT_PkRatio(ch,:) = Ch.FFT_PKS(ch,1)/Ch.FFT_PKS(ch,2);
                    Ch.wLFFT(ch,:) = w;
                    break;
                else
                    Ch.FFT_MMM(ch,:) = 0;
                    Ch.FFT_PkRatio(ch,:) = 0;
                    Ch.wLFFT(ch,:) = 0;
                end
            end
        end
        % #2 Take PSD Estimate: (Welch method)
        % Prepare static hann window: (may need to change to hann(250))
        hW = hann(windowLength);
        [Ch.PSD(ch,:), fPSD] = welch_psd(fchw(ch,:), Fs, hW);%fin-start
        % #2.2 Find Peaks and Max
        [Ch.PSD_PKS(ch,:), Ch.PSD_L(:,ch)] = findpeaks(Ch.PSD(ch,:),'SortStr','descend');
        if length(Ch.PSD_PKS(ch,:))>1
            Ch.PSD_Ltop(ch,:) = fPSD(Ch.PSD_L(1:2,ch));
            for w = 1:4
                if Ch.PSD_Ltop(ch,1)>=wLPSD(w,1) && Ch.PSD_Ltop(ch,1)<=wLPSD(w,2)
                    Ch.PSD_MMM(ch,:) = Ch.PSD_PKS(ch,1) - Ch.PSD_PKS(ch,2);
                    Ch.PSD_PkRatio(ch,:) = Ch.PSD_PKS(ch,1) / Ch.PSD_PKS(ch,2);
                    Ch.wLPSD(ch,:) = w;
                    break;
                else
                    Ch.PSD_MMM(ch,:) = 0;
                    Ch.PSD_PkRatio(ch,:) = 0;
                    Ch.wLPSD(ch,:) = 0;
                end
            end
        end
    end
    %Combine data into 'fourth' channel:
    Ch.FFT(4,:) = (Ch.FFT(1,:)+Ch.FFT(2,:)+Ch.FFT(3,:));
    [Ch.FFT_PKS(4,:), Ch.FFT_L(4,:)] = findpeaks(Ch{4}.FFT,'SortStr','descend');
    if length(Ch.FFT_PKS(4,:))>1
        Ch.FFT_Ltop(4,:) = f(Ch.FFT_L(4,1:2));
        for w = 1:4
            if Ch.FFT_Ltop(4,1)>wLFFT(w,1) && Ch.FFT_Ltop(4,1)<wLFFT(w,2) 
                Ch.FFT_MMM(4,:) = Ch.FFT_PKS(4,1) - Ch.FFT_PKS(4,2);
                Ch.FFT_PkRatio(4,:) = Ch.FFT_PKS(4,1)/Ch.FFT_PKS(4,2);
                Ch.wLFFT(4,:) = w;
                break;
            else
                Ch.FFT_MMM(4,:) = 0;
                Ch.FFT_PkRatio(4,:) = 0;
                Ch.wLFFT(4,:) = 0;
            end
        end
    end
    Ch.PSD(4,:) = Ch.PSD(1,:)+Ch.PSD(2,:)+Ch.PSD(3,:);
    [Ch.PSD_PKS(4,:), Ch.PSD_L(4,:)] = findpeaks(Ch{4}.PSD,'SortStr','descend');
    if length(Ch.PSD_PKS(4,:))>1
%         Ch.PSD_L(4,:) = Ch.PSD_L(4,:)(:);
        Ch.PSD_Ltop(4,:) = fPSD(Ch.PSD_L(4,1:2));
        for w = 1:4
            if Ch.PSD_Ltop(4,1)>=wLPSD(w,1) && Ch.PSD_Ltop(4,1)<=wLPSD(w,2)
                Ch.PSD_MMM(4,:) = Ch.PSD_PKS(4,1) - Ch.PSD_PKS(4,2);
                Ch.PSD_PkRatio(4,:) = Ch.PSD_PKS(4,1) / Ch.PSD_PKS(4,2);
                Ch.wLPSD(4,:) = w;
                break;
            else
                Ch.PSD_MMM(4,:) = 0;
                Ch.PSD_PkRatio(4,:) = 0;
                Ch.wLPSD(4,:) = 0;
            end
        end
    end
    for chn = 1:nCh+1
        FFTPeaks1(1,chn) = Ch.FFT_Ltop(chn,1);
        FFTPeaks2(1,chn) = Ch.FFT_Ltop(chn,2);
    end
    averageFFTPeak = mean([Ch.FFT_Ltop(1,1) Ch.FFT_Ltop(2,1) ...
        Ch.FFT_Ltop(3,1) Ch.FFT_Ltop(4,1)]);
    averageFFTPeak2 = mean([Ch.FFT_Ltop(1,2) Ch.FFT_Ltop(2,2) ...
        Ch.FFT_Ltop(3,2) Ch.FFT_Ltop(4,2)]);
    b1 = (Ch.wLFFT(1,:)~=0) && (Ch.wLFFT(2,:)~=0) && ...
        (Ch.wLFFT(3,:)~=0) && (Ch.wLFFT(4,:)~=0);
    averagePSDPeak = mean([Ch.PSD_Ltop(1,1) Ch.PSD_Ltop(2,1) ...
        Ch.PSD_Ltop(3,1) Ch.PSD_Ltop(4,1)]);
    b2 = (Ch.wLPSD(1,:)~=0) && (Ch.wLPSD(2,:)~=0) && (Ch.wLPSD(3,:)~=0) && (Ch.wLPSD(4,:)~=0);
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
        F(ch,:) = [Ch.FFT_Ltop(ch,1) Ch.FFT_Ltop(ch,2) Ch.FFT_MMM(ch,:) Ch.FFT_PkRatio(ch,:) ...
            Ch.wLFFT(ch,:) Ch.PSD_Ltop(ch,1) Ch.PSD_Ltop(ch,2) Ch.PSD_MMM(ch,:) ...
            Ch.PSD_PkRatio(ch,:) Ch.wLPSD(ch,:)];
    end
    Extras = [FFTPeaks1 FFTPeaks2 averageFFTPeak averageFFTPeak2 averagePSDPeak b1 b2];
    F = [F(1,:) F(2,:) F(3,:) F(4,:) Extras];
end %END FUNCTION

