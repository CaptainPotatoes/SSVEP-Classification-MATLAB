function [ Y ] = classifySSVEP2( X, plotData, thresholdFraction )
%Classify SSVEP: Y outputs likely class:
range = 250:250:length(X); % 1-4 s at 60pt intervals
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
    for d = 1:length(defFreqs)
        distFrom(d) = P(i,1) - defFreqs(d);
    end
end
[~, CLASS] = min(abs(distFrom)); 

if plotData
    commandwindow; CfLASS = input('Approve/continue?\n');
%     if isempty(CLASS)
%         if sum(b)==0
%             CLASS = 1%ClusterLoc
%         else
%             CLASS = 0
%         end
%     end
% else
%     if sum(b)==0
%         CLASS = 1
%     else
%         CLASS = 0
%     end
end
Y = CLASS;
end

