function [ CLASS ] = classifySSVEP( X, plotData, thresholdFraction )
%CLASSIFYSSVEP - FINAL VERSION FOR MATLAB CODER
    % INPUT VARS:
    % X - input array (any size)
    % start - where to start from in 'X'
    % Fs - signal sampling frequency
% range - range of window sizes to view
start = 1;
% range = 250:250:length(X); % 1-4 s at 60pt intervals
range = 1000;
NUMP = 56;
P = zeros(size(range,2),NUMP);
for i = 1:size(range,2)
    fin = start + (range(i)-1);
    fch = ssvepcfilt(X(start:fin));
    %%%Feature Extraction: (per channel)
%     fprintf('Current index = [%d to %d]\r\n',start, fin);
%     fprintf('length = %d\r\n',range(i));
%     P(i,:) = fESSVEP2(fch,false);
    [~,P(i,:)] = fESSVEP(fch,250,plotData);
end

idx = 1:4;
M = zeros(1,4);
L1 = M;
f = P(1,1:28);
if plotData
    figure(13);hold on;xlim([8 20]);
end
for i = 1:4
    [M1,L1(i)] = max(P(:,29+((i-1)*7):35+((i-1)*7)));
    M(i) = max(M1);
    if plotData
        plot(f((i-1)*7+L1(i)),M(i),'or'),xlabel('Frequency (Hz)'),ylabel('Power Density (W/Hz)'),title('Power Spectral Density Est. of Modified Signal');
    end
end
[Peak,ClusterLoc] = max(M);
idx2 = idx(ClusterLoc~=idx);
b = zeros(1,length(idx2));
Threshold = Peak/thresholdFraction;
for i=1:length(idx2)
    b(i) = (M(idx2(i)) > Threshold);
end

if plotData
    for i = 1:size(P,1)
        plot(P(i,1:28),P(i,29:end),'-*')
    end
    h = refline([0,Threshold]); h.Color = 'r';
end
% Apply other methods of classification?
if plotData
    commandwindow;CLASS = input('Approve/continue?\n');
    if isempty(CLASS)
        if sum(b)==0
            CLASS = ClusterLoc
        else
            CLASS = 0
        end
    end
else
    if sum(b)==0
        CLASS = ClusterLoc
    else
        CLASS = 0
    end
end
if plotData
    clf(13)
end
% for i = 1:size(F,1)
%     FS(1+size(F,2)*(i-1):size(F,2)*(i)) = F(i,:);
% end

end

