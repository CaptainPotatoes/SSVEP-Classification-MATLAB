function [ T, F1, S ] = extractSpectrograms( X, plotData )
%STFT Extraction:
winLim = [8,36];
Fs = 250;
wL = length(X);
nfft = 2^nextpow2(wL);
F = (0:(ceil((1+nfft)/2))-1)*Fs/nfft;
select = F<winLim(2) & F>winLim(1);
F1 = F(select);
fch = ssvepcfilt4_35(X); % filter signal
h = 32;
    wlen = 128;
if wL >= 500
    wlen = 256;  
end

K = sum(hammPeriodic(wlen))/wlen;
[Sx, ~, T] = stft2( fch, wlen, h, nfft, Fs );
S = 20*log10( abs( Sx(select,:) ) /wlen/K + 1E-6 );
if (plotData)
    figure(15); imagesc(T,F1,S),ylim(winLim),xlim([min(T),max(T)]);set(gca,'YDir','normal');colorbar;colormap(jet);
    commandwindow; a = input('Continue? \n');
end

end

