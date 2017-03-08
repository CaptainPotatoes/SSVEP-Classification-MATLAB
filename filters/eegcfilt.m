function [ Y ] = eegcfilt( X )
%EOGCFILT EEG filter for conversion to C. 
% Vectorize:
X = X(:);
% Fs = 250, N = 5
% flim = [8 18], bandpass
b = [2.13961520749732e-05,0,-0.000106980760374866,0,0.000213961520749732,0,-0.000213961520749732,0,0.000106980760374866,0,-2.13961520749732e-05];
a = [1,-8.77043379286888,35.0068378010024,-83.7229808056309,132.845833785487,-146.117834417428,112.823239428442,-60.3894491294140,21.4471017127118,-4.56451967201817,0.442209182399621];
Y = filtfilt(b,a,X);

end

