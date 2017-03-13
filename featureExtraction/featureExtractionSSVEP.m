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
wLFFT(1,:) = [9.6 10.4];%-% windows around certain target frequencies
wLFFT(2,:) = [11.9 12.7]; 
wLFFT(3,:) = [14.6 15.5];
wLFFT(4,:) = [16.2 16.7];
%----PSD----%
wLPSD(1,:) = [9.9 10.1];
wLPSD(2,:) = [12 13]; 
wLPSD(3,:) = [14.9 15.1];
wLPSD(4,:) = [16 17];
%----PREALLOCATE----%
nCh = 3;
Ch = cell(nCh+1, 1);
FFTPeaks1 = zeros(1,nCh+1);
FFTPeaks2 = zeros(1,nCh+1);
fchw = cell(nCh,1);
fchw{1,1} = fch1;
fchw{2,1} = fch2;
fchw{3,1} = fch3;
% all windows should be the same size:
windowLength = length(fch1);
% Data is already filtered:
%between 250?500dp
if windowLength>=250 && windowLength<500
    for ch=1:nCh
    % #1 Take FFT:
        [f, Ch{ch}.FFT] = get_nfft_data(fchw{ch}, Fs, 2048);
        % #1.1 Find Peaks and M/I
        [Ch{ch}.FFT_PKS, Ch{ch}.FFT_L] = findpeaks(Ch{ch}.FFT,'SortStr','descend');
        Ch{ch}.FFT_L = Ch{ch}.FFT_L(:);
        if length(Ch{ch}.FFT_PKS)>1
            %Peak max minus min
            Ch{ch}.FFT_Ltop = f(Ch{ch}.FFT_L(1:2,1));
            for w = 1:4
                if Ch{ch}.FFT_Ltop(1,1)>wLFFT(w,1) && Ch{ch}.FFT_Ltop(1,1)<wLFFT(w,2) 
                    Ch{ch}.FFT_MMM = Ch{ch}.FFT_PKS(1) - Ch{ch}.FFT_PKS(2);
                    Ch{ch}.FFT_PkRatio = Ch{ch}.FFT_PKS(1)/Ch{ch}.FFT_PKS(2);
                    Ch{ch}.wLFFT = w;
                    break;
                else
                    Ch{ch}.FFT_MMM = 0;
                    Ch{ch}.FFT_PkRatio = 0;
                    Ch{ch}.wLFFT = 0;
                end
            end
        end
        % #2 Take PSD Estimate: (Welch method)
        % Prepare static hann window: (may need to change to hann(250))
        hW = hann(windowLength);
        [Ch{ch}.PSD, fPSD] = welch_psd(fchw{ch}, Fs, hW);%fin-start
        % #2.2 Find Peaks and Max
        [Ch{ch}.PSD_PKS, Ch{ch}.PSD_L] = findpeaks(Ch{ch}.PSD,'SortStr','descend');
        Ch{ch}.PSD_L = Ch{ch}.PSD_L(:);
        if length(Ch{ch}.PSD_PKS)>1
            Ch{ch}.PSD_Ltop = fPSD(Ch{ch}.PSD_L(1:2,1));
            for w = 1:4
                if Ch{ch}.PSD_Ltop(1,1)>=wLPSD(w,1) && Ch{ch}.PSD_Ltop(1,1)<=wLPSD(w,2)
                    Ch{ch}.PSD_MMM = Ch{ch}.PSD_PKS(1) - Ch{ch}.PSD_PKS(2);
                    Ch{ch}.PSD_PkRatio = Ch{ch}.PSD_PKS(1) / Ch{ch}.PSD_PKS(2);
                    Ch{ch}.wLPSD = w;
                    break;
                else
                    Ch{ch}.PSD_MMM = 0;
                    Ch{ch}.PSD_PkRatio = 0;
                    Ch{ch}.wLPSD = 0;
                end
            end
        end
    end
    %Combine data into 'fourth' channel:
    Ch{4}.FFT = (Ch{1}.FFT+Ch{2}.FFT+Ch{3}.FFT);
    [Ch{4}.FFT_PKS, Ch{4}.FFT_L] = findpeaks(Ch{4}.FFT,'SortStr','descend');
    if length(Ch{4}.FFT_PKS)>1
        Ch{4}.FFT_L = Ch{4}.FFT_L(:);
        Ch{4}.FFT_Ltop = f(Ch{4}.FFT_L(1:2,1));
        for w = 1:4
            if Ch{4}.FFT_Ltop(1,1)>wLFFT(w,1) && Ch{4}.FFT_Ltop(1,1)<wLFFT(w,2) 
                Ch{4}.FFT_MMM = Ch{4}.FFT_PKS(1) - Ch{4}.FFT_PKS(2);
                Ch{4}.FFT_PkRatio = Ch{4}.FFT_PKS(1)/Ch{4}.FFT_PKS(2);
                Ch{4}.wLFFT = w;
                break;
            else
                Ch{4}.FFT_MMM = 0;
                Ch{4}.FFT_PkRatio = 0;
                Ch{4}.wLFFT = 0;
            end
        end
    end
    Ch{4}.PSD = Ch{1}.PSD+Ch{2}.PSD+Ch{3}.PSD;
    [Ch{4}.PSD_PKS, Ch{4}.PSD_L] = findpeaks(Ch{4}.PSD,'SortStr','descend');
    if length(Ch{4}.PSD_PKS)>1
        Ch{4}.PSD_L = Ch{4}.PSD_L(:);
        Ch{4}.PSD_Ltop = fPSD(Ch{4}.PSD_L(1:2,1));
        for w = 1:4
            if Ch{4}.PSD_Ltop(1,1)>=wLPSD(w,1) && Ch{4}.PSD_Ltop(1,1)<=wLPSD(w,2)
                Ch{4}.PSD_MMM = Ch{4}.PSD_PKS(1) - Ch{4}.PSD_PKS(2);
                Ch{4}.PSD_PkRatio = Ch{4}.PSD_PKS(1) / Ch{4}.PSD_PKS(2);
                Ch{4}.wLPSD = w;
                break;
            else
                Ch{4}.PSD_MMM = 0;
                Ch{4}.PSD_PkRatio = 0;
                Ch{4}.wLPSD = 0;
            end
        end
    end
    for chn = 1:nCh+1
        FFTPeaks1(1,chn) = Ch{chn}.FFT_Ltop(1);
        FFTPeaks2(1,chn) = Ch{chn}.FFT_Ltop(2);
    end
    averageFFTPeak = mean([Ch{1}.FFT_Ltop(1) Ch{2}.FFT_Ltop(1) ...
        Ch{3}.FFT_Ltop(1) Ch{4}.FFT_Ltop(1)]);
    averageFFTPeak2 = mean([Ch{1}.FFT_Ltop(2) Ch{2}.FFT_Ltop(2) ...
        Ch{3}.FFT_Ltop(2) Ch{4}.FFT_Ltop(2)]);
    b1 = (Ch{1}.wLFFT~=0) && (Ch{2}.wLFFT~=0) && ...
        (Ch{3}.wLFFT~=0) && (Ch{4}.wLFFT~=0);
    averagePSDPeak = mean([Ch{1}.PSD_Ltop(1) Ch{2}.PSD_Ltop(1) ...
        Ch{3}.PSD_Ltop(1) Ch{4}.PSD_Ltop(1)]);
    b2 = (Ch{1}.wLPSD~=0) && (Ch{2}.wLPSD~=0) && (Ch{3}.wLPSD~=0) && (Ch{4}.wLPSD~=0);
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
    F_1 = [Ch{1}.FFT_Ltop(1) Ch{1}.FFT_Ltop(2) Ch{1}.FFT_MMM Ch{1}.FFT_PkRatio Ch{1}.wLFFT ...
        Ch{1}.PSD_Ltop(1) Ch{1}.PSD_Ltop(2) Ch{1}.PSD_MMM Ch{1}.PSD_PkRatio Ch{1}.wLPSD]; %10 features
    F_2 = [Ch{2}.FFT_Ltop(1) Ch{2}.FFT_Ltop(2) Ch{2}.FFT_MMM Ch{2}.FFT_PkRatio Ch{2}.wLFFT ...
        Ch{2}.PSD_Ltop(1) Ch{2}.PSD_Ltop(2) Ch{2}.PSD_MMM Ch{2}.PSD_PkRatio Ch{2}.wLPSD];
    F_3 = [Ch{3}.FFT_Ltop(1) Ch{3}.FFT_Ltop(2) Ch{3}.FFT_MMM Ch{3}.FFT_PkRatio Ch{3}.wLFFT ...
        Ch{3}.PSD_Ltop(1) Ch{3}.PSD_Ltop(2) Ch{3}.PSD_MMM Ch{3}.PSD_PkRatio Ch{3}.wLPSD];
    F_4 = [Ch{4}.FFT_Ltop(1) Ch{4}.FFT_Ltop(2) Ch{4}.FFT_MMM Ch{4}.FFT_PkRatio Ch{4}.wLFFT ...
        Ch{4}.PSD_Ltop(1) Ch{4}.PSD_Ltop(2) Ch{4}.PSD_MMM Ch{4}.PSD_PkRatio Ch{4}.wLPSD];
    Extras = [FFTPeaks1 FFTPeaks2 averageFFTPeak averageFFTPeak2 averagePSDPeak b1 b2];
    F = [F_1 F_2 F_3 F_4 Extras];
end %END FUNCTION

