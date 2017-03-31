function [ TP, OUT ] = featureAnalysis(FEA0, len)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% TP = TOTAL POINTS
threshold = 5.5;
FREQ = [10 12 15 16];
TP = zeros(1,4);
OUT = 0;
O10 = FEA0(1:4,:);
O12 = FEA0(5:8,:);
O15 = FEA0(9:12,:);
O16 = FEA0(13:16,:);
S1 = [sum(O10(:)==0) sum(O12(:)==0) sum(O15(:)==0) sum(O16(:)==0)];
P1 = 24-S1;
% PSD
%                     Extract PSD STUFF:
PSD10 = FEA0(1:4,3:4);
PSD12 = FEA0(5:8,3:4);
PSD15 = FEA0(9:12,3:4);
PSD16 = FEA0(13:16,3:4);
S2 = [sum(PSD10(:)==0) sum(PSD12(:)==0) sum(PSD15(:)==0) sum(PSD16(:)==0)];
P2 = 8-S2;
varFreqs = [var(PSD10(:,1)) var(PSD12(:,1)) var(PSD15(:,1)) var(PSD16(:,1)) ]
AvgMax = [mean(PSD10(:,2)) mean(PSD12(:,2)) mean(PSD15(:,2)) mean(PSD16(:,2))]
minAvgMax = min(AvgMax(AvgMax~=0));
if ~isempty(minAvgMax)
    compare = AvgMax./minAvgMax
else
    compare = zeros(1,4)
end
% TREE:
for i=1:4
    if(P2(i)~=0)
        TP(1,i) = compare(i); %TODO
    end
end
TP 
SortTP = sort(TP,'descend');
[M, I] = max(TP);
Dif = M - SortTP(2) %Difference between first and second place
% TODO: CONDITIONS FOR TRACKING HISTORY?
threshold = 1;
if Dif > threshold
    % select possible signal.
    if varFreqs(I) < 0.1
        OUT = FREQ(I)
    end
end
% if varFreqs(I) < 0.1
%     OUT = FREQ(I)
% end
end

