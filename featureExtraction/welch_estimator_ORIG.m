%% Function for spectra estimation by Welch's method
% Developed by Luiz A. Baccala, Fl?vio Caduda and Luciano Caldas, all from
% Escola Polit?cnica - Poli-USP, with cooperation of Carlos Pagani and Felipe
% Amaral from Escola de Engenharia de S?o Carlos - EESC-USP.
%
% Cross-spectra matrix are estimated by Welch's method with 50% overlap and
% the window energy loss are compasated by a factor of 1/sum(Wi.^2) where
% Wi are the elements of the window [1]. Then, the spectra becomes:
% Sxy = fft(x)*conj(fft(y))/sum(Wi.^2)
%
% Code was tested with a known- spectra signal from a white noise filtered
% by a filter. The variance (power) of the signal checks with the integral
% of the PSD estimated.
%
% INPUT:
% -- signals: matrix of signals to perform the spectra estimatino. Size is
% [Samples x number of sensors];
% -- fs: samplerate in Hertz;
% -- window: data taper desired. Must be a vector. For best performance it
% should be a power of 2. For general applications do: window=hanning(1024);
%
% OUTPUT:
% -- CSM: Cross Spectral Matrix: Unilateral (0:fs/2) spectra. Welch's
% method is used with 50% overlap. Matrix size: sensors x sensors x
% windowsize/2
% -- frequencies: vector with all frequencies corresponding to each layer
% (3rd layer in depth) of CSM.
%
% LAST REVISION: Aug - 18 - 2016
% ADDED 'fs' missing term in line 82, for calibration factor
% [1] Trobs,M.; Heinzel,G. "Improved spectrum estimation from digitized
% time series on a logarithmic frequency axis"
% doi:10.1016/j.measurement.2005.10.010
function [CSM,frequencies] = welch_estimator(signals,fs,window)
if size(signals,2) > size(signals,1)
    signals = signals.';
end
window = window(:);
sensors = size(signals,2);
windowsize = length(window);
frequencies = (0:(windowsize/2-1))*fs/windowsize;
block_samples = windowsize; %must be even, best if 2^n
signal_samples = size(signals,1);
number_of_signals = size(signals,2);
back_shift = block_samples./2; %ORIGINAL;
number_of_blocks = floor((2*signal_samples)./block_samples) -1;
data_taper = window;
data_taper = repmat(data_taper,1,number_of_signals);
% Data segmentation into blocks of size block_samples:
S = zeros(block_samples/2,number_of_signals.^2); %ORIGINAL
% S = zeros(ceil(block_samples/2),number_of_signals.^2);
for a = 1:number_of_blocks
    % Retrieve current data block
    Data_Block = signals((a-1)*back_shift+1:block_samples +(a-1)*back_shift,:);
    Data_Block = Data_Block - repmat(mean(Data_Block),block_samples,1);
    Data_Block = Data_Block.*data_taper; %Taper it
    Data_Block = fft(Data_Block); %FFT it,
    % bilateral DFT
    % viii
    Data_Block = Data_Block(1:block_samples/2,:); %ORIGINAL
    % Data_Block = Data_Block(1:ceil(block_samples/2),:);
    %All spectral combinations:
    P = zeros(block_samples/2,number_of_signals.^2); %ORIGINAL
    % P = zeros(ceil(block_samples/2)/2,number_of_signals.^2);
    c = 1;
    for aa = 1:size(Data_Block,2)
        for b = aa:size(Data_Block,2)
            % P(:,c) = Data_Block(:,b).*conj(Data_Block(:,aa)); % THIS
            % IS FOR WIND TUNNEL EESC-USP BEAMFORMING CODE
            % P(:,c) = Data_Block(:,aa).*conj(Data_Block(:,b)); % THIS IS THE ORIGINAL
            % LINE
            P(:,c) = real(Data_Block(:,aa).*conj(Data_Block(:,b)));
            % IS FOR FAN RIG BEAMFORMING CODE
            c = c+1;
        end
    end
    % Sum the spectrums up ...
    S = S + P;
end
S = S*2/(sum(window.^2)*fs*number_of_blocks); % Average them out
Sf = zeros(sensors,sensors,size(S,1));
c=1;
for a = 1:sensors
    for b = a:sensors
        Sf(a,b,:) = S(:,c);
        c = c+1;
    end
end
% clear S
CSM = Sf;
for i = 1:size(CSM,3)
    CSM(:,:,i) = CSM(:,:,i) + CSM(:,:,i)' - eye(sensors).*CSM(:,:,i);
end

end