function [ PSD ] = extractConvolutedPowerSpectrum( X1, X2, plotData)

fch1 = ssvepcfilt2(X1(1:end)); %[5 40]
fch2 = ssvepcfilt2(X2(1:end));
conv2ch = conv(fch1,fch2,'full');

X = conv2ch(1:end-1);
X = X(:);

wL = length(X);
hW = hannWin(wL);

[PSD, fPSD] = welch_psd(X, 250, hW);

if plotData
    figure(13);hold on;xlim([8 30]);
    plot(fPSD, PSD);
end

if plotData
%     commandwindow; IRINPUT = input('Approve/continue?\n');
end
if plotData
    clf(13)
end
end



