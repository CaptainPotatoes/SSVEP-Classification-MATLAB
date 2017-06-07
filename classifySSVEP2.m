function [ CLASS, F, M ] = classifySSVEP2( X, plotData, thresholdFraction )
%CLASSIFYSSVEP - FINAL VERSION FOR MATLAB CODER - thresholdFraction 
    % INPUT VARS:
    % X - input array (any size)
    % start - where to start from in 'X'
    % Fs - signal sampling frequency
% range - range of window sizes to view
start = 1;
% range = 250:250:length(X); % 1-4 s at 60pt intervals
range = 1000;
NUMP = 56;
% P = zeros(size(range,2),NUMP);
for i = 1:size(range,2)
    fin = start + (range(i)-1);
    fch = ssvepcfilt2(X(start:fin)); %[5 40]
%     fprintf('Current index = [%d to %d]\r\n',start, fin);
%     fprintf('length = %d\r\n',range(i));
%     [F, P(i,:), f(i,:)] = fESSVEP(fch,250,plotData); % Extract Features
    
end
ClusterSize = length(P)/4;
idx = 1:4;
% M = zeros(1,4);
if plotData
    figure(13);hold on;xlim([8 20]);
%     for i=1:length(range)
%         plot(f(i,:),P(i,:));
%     end
end
for j = 1:length(range)
    plot(f(j,:),P(j,:));
    for i = 1:4
        start = (i-1)*ClusterSize + 1;
        fin = start+(ClusterSize-1); 
        [M1,L1] = max(P(j,start:fin));
        [M(j,i)] = max(M1);
        if plotData
            plot(f(j,(i-1)*ClusterSize+L1),M(j,i),'or'),xlabel('Frequency (Hz)'),ylabel('Power Density (W/Hz)'),title('Power Spectral Density Est. of Modified Signal');
        end
    end
    [Peak(j),ClusterLoc(j)] = max(M(j,:));
    Threshold(j) = Peak(j)/thresholdFraction;
    idx2(j,:) = idx(ClusterLoc(j)~=idx);
    for k=1:size(idx2,2)
        b(j,k) = (M(j,idx2(k)) > Threshold(j));
    end
end

if plotData
    h = refline([0,max(Threshold)]); h.Color = 'r';
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

