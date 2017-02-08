function [ Y ] = ssvepFilter( X, Fs, Fpoint, Frange, NOrder )
% Modular SSVEP Filter 
% The idea is to pass (data, sampling freq, and the frequency you are
% looking for, and this function will narrowly filter around it. 

% X is the data window
% Fs is sampling frequency
% Fpoint is the sampling point of inspection. 
% Frange is half the size of filter window (e.g. +/- Frange)
% NOrder is the order of the filter

% Y (returned) is the filtered data. 

if ~isvector(X)
    error('Input data not a vector');
end

if Fpoint > 4
    f1 = Fpoint - Frange;
    f2 = Fpoint + Frange;
    Wn = [f1 f2]*2/Fs;
    N=NOrder;
    [a,b] = butter(N,Wn,'bandpass'); %bandpass filtering
    Y = filtfilt(a,b,X);
else 
    error('Fpoint must be a real number greater than 4');
end

end

