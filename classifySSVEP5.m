function [ Y , CLASS0 ] = classifySSVEP5( X1, X2, plotData, thresholdFraction )
start = 1; % start = 501;
CLASS_LABELS = [1,2,3,4];
range = length(X1);
for i = 1:size(range,2)
    fch1 = ssvepcfilt2(X1(start:end)); %[5 40]
    fch2 = ssvepcfilt2(X2(start:end));
    conv2ch = conv(fch1,fch2,'full');
    [Ppsd] = fPowerSpectrum5(conv2ch(1:end-1),250,plotData);
end
if plotData
    figure(13);hold on;xlim([8 30]);
end

[M, Idx] = max(Ppsd);
%% Version 1: (Simple)
CLASS0 = CLASS_LABELS(Idx)
%% Version 2: (More complex, uses threshold):
% %{
IF = 1:4;
Idx2 = IF((1:4)~=Idx);
b = zeros(1,length(Idx2));
Threshold = M/thresholdFraction;
for i=1:length(Idx2)
    b(i) = (Ppsd(Idx2(i)) > Threshold);
end
if sum(b)==0
    CLASS = CLASS_LABELS(Idx);
else
    CLASS = 0;
end
%}

if plotData
    h = refline([0,Threshold]); h.Color = 'r'; commandwindow; IRINPUT = input('Approve/continue?\n');
end
Y = CLASS;
if plotData
    clf(13)
end
end

