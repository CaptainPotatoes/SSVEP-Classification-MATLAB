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
%%%%%%---- COPY FROM 'TEMP' WHEN DONE, REMOVE ALL PLOTS AND FPRINTFs
%%%%%% --- TEST WITH functionTestEOG->fullHybridClassifier
end %END FUNCTION

