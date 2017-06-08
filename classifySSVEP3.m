function [ Y ] = classifySSVEP3( X1, X2, plotData, thresholdFraction )
start = 1;
% range = 500:250:length(X); % 1-4 s at 60pt intervals
range = 1000;
NUMP = 56;
P = zeros(size(range,2),NUMP);
for i = 1:size(range,2)
    fin = start + (range(i)-1);
    fch = ssvepcfilt2(X1(start:fin)); %[5 40]
    fch2 = ssvepcfilt2(X2(start:fin));
    conv2ch = conv(fch,fch2,'full'); 
%     if  mod(length(conv2ch),2)==1
        P(i,:) = fECONV2(conv2ch(1:end-1),250,plotData);
%     else
%         P(i,:) = fECONV2(conv2ch,250,plotData);
%     end
end

idx = 1:4;
M = zeros(1,4);
L1 = M;
f = P(1,1:28);
if plotData
    figure(13);hold on;xlim([8 30]);
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
if (nargin<5)
    Threshold = Peak/thresholdFraction;
end
%TODO: 7.3E-13
for i=1:length(idx2)
    b(i) = (M(idx2(i)) > Threshold);
end

if plotData
%     type = ['.r','.b','.k','.c']; 
    for i = 1:size(P,1)
%         plot(P(i,1:28),P(i,29:end),'.r');
        for j = 1:4
            plot(P(i,((j-1)*7 + 1):j*7),P(i,((j+3)*7 + 1):(j+4)*7),'-.');
        end
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
Y = CLASS;
if plotData
    clf(13)
end
end
