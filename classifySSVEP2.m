function [ Y ] = classifySSVEP2( X, plotData )
%Classify SSVEP: Y outputs likely class:
range = 1000;
start = 1;
defFreqs = [10, 12.50, 15.15, 16.66];
% range = 1000;
% NUMP = 56;
% P = zeros(size(range,2),NUMP);
for i = 1:size(range,2)
    fin = start + (range(i)-1);
    fch = ssvepcfilt(X(start:fin));
    %%%Feature Extraction: (per channel)
%     fprintf('Current index = [%d to %d]\r\n',start, fin);
%     fprintf('length = %d\r\n',range(i));
    P(i,:) = fESSVEP3(fch,plotData);
%     for d = 1:length(defFreqs)
%         distFrom(d) = P(i,1) - defFreqs(d);
%     end
end

if plotData
    commandwindow; CfLASS = input('Approve/continue?\n');
end
Y = 0;
end

