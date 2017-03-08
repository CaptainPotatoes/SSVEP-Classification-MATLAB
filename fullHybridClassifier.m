function [ Y ] = fullHybridClassifier( ch1, ch2, ch3, ch4, tXEOG, tYEOG)
%Full hybrid EOG/EEG classifier
% Ch1 = Fp1
% Ch2 = Fp2
% Ch3 = Fpz
% Ch4 = Right Eye
% 
% tXEOG = Training Data for EOG
% tYEOG = Classes for EOG
% 
% tX    = SSVEP Training Data (not sure if I will include just yet).
% tY    = SSVEP Training Classes
% 
% The first part of this classifier determines if an emergency stop is
% requested in the form of a double-blink by the subject:
% A '1' is a double blink.
% Any other result is a pass, and classification continues.
% The SSVEP Portion will output one of the following corresponding to its
% frequency makeup:
% CLASS :: CORRESPONDING FREQ (ACTUAL)
% 10    :: 10.00Hz
% 12    :: 12.50Hz
% 15    :: 15.15Hz
% 16    :: 16.66Hz

% Window length??
Y=0; %Default Value.
chLen = length(ch1);
numFeatures = 10;
fch1f = zeros(size(ch1,2),numFeatures);
fch2f = zeros(size(ch1,2),numFeatures);
fch3f = zeros(size(ch1,2),numFeatures);
fch4f = zeros(size(ch1,2),numFeatures);
DB = false;
SSVEP_PRESENT = false;
if chLen>=250
    for i = 1:size(ch1,2);
        % Filter using optimized EOG filter: 
        % (take last second, regardless of actual window length):
        ch1f = eogcfilt(ch1(end-249:end,i));
        ch2f = eogcfilt(ch2(end-249:end,i));
        ch3f = eogcfilt(ch3(end-249:end,i));
        ch4f = eogcfilt(ch4(end-249:end,i));
        %Extract EOG features: (1s window)
        fch1f(i,:) = featureExtractionEOG(ch1f');
        fch2f(i,:) = featureExtractionEOG(ch2f');
        fch3f(i,:) = featureExtractionEOG(ch3f');
        fch4f(i,:) = featureExtractionEOG(ch4f');
    end
    %Combine features:
    samplesX = [fch1f fch2f fch3f fch4f] ;
    %Boolean DB: represents presence of a double blink.
    Y = knn(samplesX, tXEOG, tYEOG, 3);
    if Y==1
        DB = true;
    end
end
% SSVEP CLASSIFICATION: 
% PRECONDITIONS: EOG must not have been triggered. 
    % starts with 1/2 second analysis and moves up. 
    % Output can be one of the four SSVEP classes [10 12 15 16]

if chLen>=124
    if ~DB
        %If no double blink has been detected in final second of data. 
        % Use a decision tree.
        %Y = knn(tsX, tX, tY, 1); %Fine KNN
        ch1f = eegcfilt(ch1);
        ch2f = eegcfilt(ch2);
        ch3f = eegcfilt(ch3);
        ch4f = eegcfilt(ch4);
        % Extract EEG Features
        SSVEPfch1(i,:) = featureExtractionSSVEP(ch1f');
        SSVEPfch2(i,:) = featureExtractionSSVEP(ch2f');
        SSVEPfch3(i,:) = featureExtractionSSVEP(ch3f');
        SSVEPfch4(i,:) = featureExtractionSSVEP(ch4f');
        % Analysis:
        samplesX = [SSVEPfch1 SSVEPfch2 SSVEPfch3 SSVEPfch4] ;
        % Decisions
        if ~SSVEP_PRESENT
            Y = 0;
        end
    end
end


end

