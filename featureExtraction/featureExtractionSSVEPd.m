function [ f, FFT, fPSD, PSD ] = featureExtractionSSVEPd( X, Fs )
%featureExtraction - DEPRECATED!
% For single channel feature extraction. 
% X - new samples (SINGLE CHANNEL; FILTERED IN FULL CLASSIFIER)

% Preload thresholds:
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
% Check length
windowLength = length(X);
%between 250?500dp
if windowLength>=250 && windowLength<500
    %Classification method #1 (w/o STFT):
    % #1.0 Take FFT:
    [f, FFT] = get_nfft_data(X(1:250), Fs, 2048);
    % #1.1 Find any peaks:
    [FFT_PKS, FFT_L] = findpeaks(FFT,'SortStr','descend');
    % #1.2 Convert to column vector
    FFT_L = FFT_L(:);
    % #1.3 If there are more than one peaks:
    if(length(FFT_PKS)>1)
        % Select the top two corresponding frequencies:
        FFT_Ltop = f(FFT_L(1:2,1));
        % iterate through thresholds (wLFFT) and, if FFT_Ltop(1) is between
        % them, save max-min data as well as max/min ratio data. Saves 'w'
        % index
        for w = 1:4
            if FFT_Ltop(1,1)>wLFFT(w,1) && FFT_Ltop(1,1)<wLFFT(w,2)
                FFT_MMM = FFT_PKS(1) - FFT_PKS(2);
                FFT_PkRatio = FFT_PKS(1)/FFT_PKS(2);
                wLFFT = w;
                break;
            else
                FFT_MMM = 0;
                FFT_PkRatio = 0;
                wLFFT = 0;
            end
        end
    end
    
    % #2.0 Take PSD Estimate using Welch's Method:
    window = hann(250);
    [PSD, fPSD] = welch_psd(X(1:250), Fs, window);
    % #2.1 Find Peaks:
    [PSD_PKS, PSD_L] = findpeaks(PSD,'SortStr','descend');
    % #2.2 Convert to column vector:
    PSD_L = PSD_L(:);
    if length(PSD_PKS)>1
        PSD_Ltop = fPSD(PSD_L(1:2));
        for w = 1:4
            if PSD_Ltop(1,1)>=wLPSD(w,1) && PSD_Ltop(1,1)<=wLPSD(w,2)
                PSD_MMM = PSD_PKS(1) - PSD_PKS(2);
                PSD_PkRatio = PSD_PKS(1)/PSD_PKS(2);
                wLPSD = w;
            else
                PSD_MMM = 0;
                PSD_PkRatio = 0;
                wLPSD = 0;
            end
        end
    end
elseif windowLength>=500
    %Classification method #2 (w/ STFT):
else
    error('Not enough data!\n');
end
% F = horzcat();

end