function [ Ppsd ] = fPowerSpectrum( X0, Fs, plotData )
%FESSVEP Feature Extraction for single (m x 1) SSVEP EEG Data Vector
%   X (m x 1) vectorize input:
% Fix X size:
X = zeros(1,length(X0));
X = (X0(:)');
% Fs        = Sampling Frequency;
% plotData  = Plot Data

%%%%% - Thresholds: - %%%%%
NUMBER_CLASSES = 5;
threshPSD = zeros(NUMBER_CLASSES,2);
threshPSD(1,:) = [9 13];
threshPSD(2,:) = [14 15.5];
threshPSD(3,:) = [16.1 17.2];
threshPSD(4,:) = [18.0 18.8];
threshPSD(5,:) = [19.4 20.4];
%%% - Constants - %%%
selc = ['.m';'.b';'.m';'.k';'.c']; %select dot color; 
winLim = [8,24];
% - Variables - %
if plotData
    fH = figure(12); %-% Figure Handle
    set(fH, 'Position', [0, 0, 1440, 960]);
    clf(fH)
end
wL = length(X);
if mod(wL,2) == 1
    PSD = zeros(1,(wL-1)/2);
else
    PSD = zeros(1,wL/2);
end

Lpsd = zeros(1,NUMBER_CLASSES);
Ppsd = zeros(1,NUMBER_CLASSES);

hW = hannWin(wL);
if wL >= 250
    [PSD, fPSD] = welch_psd(X, Fs, hW);
    if plotData
        hold on;plot(fPSD,PSD),xlim(winLim);
    end
    for i=1:NUMBER_CLASSES
        [fselect2, psdselect, Lpsd(i), Ppsd(i)] = get_psd_features(fPSD,PSD,threshPSD(i,:));
        if plotData
            hold on;plot(fselect2,psdselect, selc(i,:)); plot(Lpsd(i),Ppsd(i),'or');
        end
    end
end
% if plotData
%     figure(13); hold on; plot(Lpsd, Ppsd, '*');
% end


end

