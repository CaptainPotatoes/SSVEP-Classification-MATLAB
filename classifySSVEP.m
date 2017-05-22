function [ CLASS ] = classifySSVEP( X, plotData )
%CLASSIFYSSVEP - FINAL VERSION FOR MATLAB CODER
    % INPUT VARS:
    % X - input array (any size)
    % start - where to start from in 'X'
    % Fs - signal sampling frequency
    
% range - range of window sizes to view
start = 1;
range = 250:250:1000; % 1-4 s at 60pt intervals
NUMP = 56;
P = zeros(size(range,2),NUMP);
for i = 1:size(range,2)
    fin = start + (range(i)-1);
    fch = ssvepcfilt(X(start:fin));
    %%%Feature Extraction: (per channel)
%     fprintf('Current index = [%d to %d]\r\n',start, fin);
%     fprintf('length = %d\r\n',range(i));
    P(i,:) = fESSVEP2(fch,false);
end

idx = 1:4;
M = zeros(1,4);
L = M;
f = P(1,1:28);
if plotData
    figure(13);hold on;xlim([8 20]);
end
for i = 1:4
    [M(i),L(i)] = max(max(P(:,29+((i-1)*7):35+((i-1)*7))));
    if plotData
        plot(f((i-1)*7+L(i)),M(i),'or');
    end
end
[Peak,ClusterLoc] = max(M);
idx2 = idx(ClusterLoc~=idx);
b = zeros(1,length(idx2));
Threshold = Peak/2;
for i=1:length(idx2)
    b(i) = (M(idx2(i)) > Threshold);
end

if plotData
    for i = 1:size(P,1)
        plot(P(i,1:28),P(i,29:end),'*k')
    end
    h = refline([0,Threshold]); h.Color = 'r';
end

% commandwindow;CLASS = input('Approve/continue?\n');
% if isempty(CLASS)
    if sum(b)==0
        CLASS = ClusterLoc;
    else
        CLASS = 0;
    end
% end
if plotData
    clf(13)
end
% for i = 1:size(F,1)
%     FS(1+size(F,2)*(i-1):size(F,2)*(i)) = F(i,:);
% end

end

