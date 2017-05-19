function [ F ] = featureExtractionSSVEP( X, range, filtRange, cont, start, Fs )
%featureExtractionSSVEP
    % Input Vars:
%  X            = Signal Vector;
%  range        = size of windows to test (e.g. 250 ? 2500 in 60 point increments)
%  filtRange    = Bandpass range for butterworth filter
%  cont         = Plot Graph if Empty, else run to completion
%  start        = Data point to start analyzing dataset from
%  Fs           = Sampling Rate (250 Hz by default)

    % Default Values:
if nargin < 6
    Fs = 250;
end
if nargin < 5
    start = 1;
end
if nargin < 4
    cont = 0;
end
if nargin < 3
    filtRange = [8,20];
end
if nargin < 2
    range = 250:60:2500;
end

    %Preallocate:
% NUMBERFEATURES = 24;
% F = zeros(size(range,2),NUMBERFEATURES);
for i = 1:size(range,2)
    fin = start + (range(i)-1);
    fprintf('Current index = [%d to %d]\r\n',start, fin);
    fprintf('length = %d\r\n',range(i));
    fch = customFilt(X(start:fin),Fs,filtRange,3);
        %%%Feature Extraction: (per channel)
    F(i,:) = fESSVEP(fch,Fs,isempty(cont));
    if isempty(cont)
        commandwindow;
        cont = input('Approve/continue?\n');
        clf(12);
    end
end

end

