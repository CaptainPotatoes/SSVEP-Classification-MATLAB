function [ b, a ] = customFiltDesign( Fs,f,N )
%FILT_CUSTOM Allows for quick customization of bandpass filter parameters

Wn=[f(1) f(2)]*2/Fs; % cut off based on Fs
[b,a] = butter(N,Wn);

end

