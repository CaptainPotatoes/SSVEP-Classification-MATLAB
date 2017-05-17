function [ Y, T ] = testSignal( freq, len, amplitude, Fs )
%testSignal For generating idealized SSVEP Signals (Fourier Series)
%   freq = frequency of signal waveform
%   Fs = Sampling Frequency

if(nargin<3)
    amplitude = 1E-4;
end

if nargin < 4 %standard sampling freq
    Fs = 250;
end

h = 1/Fs;
Tend = len/Fs-h;
T = 0:h:Tend; %Time Signal With Specified Frequency

Y = amplitude*sin(2*pi*freq*T);
Y = Y(:);
end

