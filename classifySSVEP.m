function [ Y, Ppsd, featurePoint ] = classifySSVEP( X1, X2, plotData, thresholdFraction, start )
CLASS_LABELS = [0,1,2,3,4];
range = length(X1);
Y = zeros(2,1);
for i = 1:size(range,2)
    fch1 = ssvepcfilt2(X1(start:end)); %[5 40]
    fch2 = ssvepcfilt2(X2(start:end));
    conv2ch = conv(fch1,fch2,'full');
    if plotData
%         fH = figure(2);
%         set(fH, 'Position', [0, 0, 650*2, 540/2]);
%         clf(fH);
%         tX = 0:1/250:(length(X1)/250-1/250);
%         subplot(1,2,1); plot(tX,X1); title('Raw Data Channel 1'); xlabel('Time (s)'); ylabel('Signal Amplitude (V)');
%         subplot(1,2,2); plot(tX,X2); title('Raw Data Channel 2'); xlabel('Time (s)'); ylabel('Signal Amplitude (V)');
    end
%     [Ppsd] = fESSVEP(fch1,250,plotData,10);
%     [Ppsd] = fESSVEP(fch2,250,plotData,11);
    [Ppsd, Lpsd] = fESSVEP(conv2ch(1:end-1),250,plotData);
end
if plotData
    figure(13);hold on;xlim([8 30]);
end

[M, Idx] = max(Ppsd);
featurePoint = [Lpsd(Idx), M];
%% Version 1: (Simple)
CLASS0 = CLASS_LABELS(Idx)
%% Version 2: (More complex, uses threshold):
% %{
% if(Idx == 1) %if 0, then check that it towers over any SSVEP signals.  
%     A=1:5;
%     IDX2 = A((1:5)~=Idx);
%     M2 = M/2;
%     for i=1:length(Idx2)
%         b(i) = (Ppsd(Idx2(i)) > Threshold);
%     end
%     if sum(b)==0
%         CLASS = CLASS_LABELS(Idx)
%     else
%         CLASS = 0
%     end
% end
IF = 1:5;
Idx2 = IF((1:5)~=Idx);
b = zeros(1,length(Idx2));
Threshold = M/thresholdFraction;
for i=1:length(Idx2)
    b(i) = (Ppsd(Idx2(i)) > Threshold);
end
if sum(b)==0
    CLASS = CLASS_LABELS(Idx)
else
    CLASS = 0
end
%}

if plotData
%     h = refline([0,Threshold]); h.Color = 'r'; commandwindow; IRINPUT = input('Approve/continue?\n');
end
Y(1) = CLASS0;
Y(2) = CLASS;
if plotData
    clf(13)
end
end

