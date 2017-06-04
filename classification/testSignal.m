function [ Y, T, Y1, Y2 ] = testSignal( freq, len, amplitude, Fs )
%testSignal For generating idealized SSVEP Signals (Fourier Series)
%   freq = frequency of signal waveform
%   Fs = Sampling Frequency
%%% Outputs:
% Y = Fundamental Frequency
% T = time singal
% Y1 = first harmonic
% Y2 = second harmonic
% SQW = square wave

if(nargin<3)
    amplitude = 1E-4;
end
if nargin < 4 %standard sampling freq
    Fs = 250;
end

h = 1/Fs;
Tend = len/Fs-h;
T = 0:h:Tend; %Time Signal With Specified Frequency

Y  = amplitude*sin(2*pi*freq*T) + amplitude*cos(2*pi*freq*T);
Y  = Y(:);
Y1 = amplitude*sin(4*pi*freq*T);
Y1 = Y1(:);
Y2 = amplitude*sin(6*pi*freq*T);
Y2 = Y2(:);
% SQW = amplitude*square(2*pi*freq*T);
% SQW = SQW(:);
end

