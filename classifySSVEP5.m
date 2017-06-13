function [ Y , P_data ] = classifySSVEP5( X1, X2, plotData, thresholdFraction )
start = 1;
CLASS_LABELS = [0,1,2,3,4];
range = length(X1);
NUMP = 70;
P = zeros(size(range,2),NUMP);
for i = 1:size(range,2)
    fin = start + (range(i)-1);
    fch = ssvepcfilt2(X1(start:fin)); %[5 40]
    fch2 = ssvepcfilt2(X2(start:fin));
    conv2ch = conv(fch,fch2,'full');
%     [P, C] = fESSVEP(conv2ch(1:end-1),250,plotData);
    P = fESSVEP2(conv2ch(1:end-1),250,plotData);
end
P_len = size(P,2);
NumClasses = (P_len/2/7);
idx = 1:NumClasses;
M = zeros(1,NumClasses);
L1 = M;
f = P(1,1:(P_len/2));
P_data = P(1,(P_len/2+1):(P_len));
if plotData
    figure(13);hold on;xlim([8 30]);
end
for i = 1:NumClasses
    [M1, L1(i)] = max(P_data(:,(i-1)*7+1:7*i));
    M(i) = max(M1);
    if plotData
        plot(f((i-1)*7+L1(i)),M(i),'or'),xlabel('Frequency (Hz)'),ylabel('Power Density (W/Hz)'),title('Power Spectral Density Est. of Modified Signal');
    end
end

[Peak,ClusterLoc] = max(M);
idx2 = idx(ClusterLoc~=idx);
b = zeros(1,length(idx2));
if (nargin<5)
    Threshold = Peak/thresholdFraction;
end

for i=1:length(idx2)
    b(i) = (M(idx2(i)) > Threshold);
end

if plotData
    for i = 1:NumClasses
        plot(f((i-1)*7+1:7*i),P_data((i-1)*7+1:7*i),'*');
    end
    h = refline([0,Threshold]); h.Color = 'r';
end
% Apply other methods of classification?
if plotData
    commandwindow; CLASS = input('Approve/continue?\n');
    if isempty(CLASS)
        if sum(b)==0
            CLASS = CLASS_LABELS(ClusterLoc)
        else
            CLASS = 0
        end
    end
else
    if sum(b)==0
        CLASS = CLASS_LABELS(ClusterLoc)
    else
        CLASS = 0
    end
end
Y = CLASS;
if plotData
    clf(13)
end
end

