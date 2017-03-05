function [ Y ] = fullHybridClassifier( ch1, ch2, ch3, ch4, tX, tY )
%Full hybrid EOG/EEG classifier
% Ch 1 = Fp1
% Ch 2 = Fp2
% Ch 3 = Fpz
% Ch 4 = Right Eye
% The first part of this classifier determines if an emergency stop is
% requested in the form of a double-blink by the subject:
% A '1' is a double blink.
% Any other result is a pass, and classification continues.
% Window length??

chLen = length(ch1);
if chLen>=250
    % Make sure all channels are arrays:
    ch1 = ch1(:);
    ch2 = ch2(:);
    ch3 = ch3(:);
    ch4 = ch4(:);
    % Filter using optimized EOG filter: 
    % (take last second, regardless of actual window length):
    ch1f = eogcfilt(ch1(end-249:end));
    ch2f = eogcfilt(ch2(end-249:end));
    ch3f = eogcfilt(ch3(end-249:end));
    ch4f = eogcfilt(ch4(end-249:end));
    %Extract EOG features: (1s window)
    fch1f = featureExtractionEOG(ch1f);
    fch2f = featureExtractionEOG(ch2f);
    fch3f = featureExtractionEOG(ch3f);
    fch4f = featureExtractionEOG(ch4f);
    %Combine features:
    samplesX = [fch1f fch2f fch3f fch4f] ;
    %Boolean DB: represents presence of a double blink.
    DB = false;
    Y = knn(samplesX, tX, tY, 3);
    if Y==1
        DB = true
    end
end
% SSVEP CLASSIFICATION: 
% PRECONDITIONS: EOG must not have been triggered.
    % starts with 1/2 second analysis and moves up. 
if ~DB
    % if Y~=1
end

end

