function [ Y , F ] = fullHybridClassifier( ch1, ch2, ch3, ch4, tXEOG, tYEOG, Fs )
%Full hybrid EOG/EEG classifier
% Ch1 = Fp1
% Ch2 = Fp2
% Ch3 = Fpz
% Ch4 = Right Eye
% 
% tXEOG = Training Data for EOG
% tYEOG = Classes for EOG
% 
% tXSSVEP    = SSVEP Training Data (not sure if I will include just yet).
% tYSSVEP    = SSVEP Training Classes
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
%Temporary:
F = zeros(1,53);
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

if chLen>=250
    if ~DB
        %If no double blink has been detected in final second of data. 
        % Use a decision tree.
        %Y = knn(tsX, tX, tY, 1); %Fine KNN
        %limit size:
        ln = min([length(ch1) length(ch2) length(ch3)]);
        ch1 = ch1(1:ln);
        ch2 = ch2(1:ln);
        ch3 = ch3(1:ln);
        %Filter:
        ch1f = eegcfilt(ch1);
        ch2f = eegcfilt(ch2);
        ch3f = eegcfilt(ch3);
        % Extract SSVEP Features (Part 1 from individual channels):
        plotData = false;
        F = featureExtractionSSVEP(ch1f, ch2f, ch3f, Fs);
%         numFeats = size(F,1);
        %OR USE chLen<500; chLen>=500
        if chLen<500%numFeats <= 53
            % window length is less than 500 samples
%             Y = knn(F,tXSSVEP,tYSSVEP,5);
        elseif chLen>=500%numFeats>53
            % Window length is 500 or more samples.
%             Y = knn(F,tXSSVEPlong,tYSSVEPlong,5);
        end
        % Analysis: Use Tree-based classification or CCA-FKNN:
        %%% ^ TODO: use self-programmed decision tree for l = 250;
        %%% ^ Use CCA-FKNN for larger datasets (for l>=500);
%         numFeatures = length(F); 
        % Decisions
        if ~SSVEP_PRESENT %TEMP (obviously)
            Y = 0;
        end
    end
end


end

