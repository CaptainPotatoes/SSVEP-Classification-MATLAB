function [ Y ] = treeClassifier( F, chLen )
%treeClassifier - For SSVEP Classification. 
%  
% 
CLASS = [10 12 15 16];
Y = zeros(5,1);
if length(F)==30 %short classifier
    %re-establish vars:
    wLFFT = F(1:4);
    wLPSD = F(5:8);
    FFT_PkRatio = F(9:12);
    averagePkRatioFFT = mean(FFT_PkRatio);
    PSD_PkRatio = F(13:16);
    averagePkRatioPSD = mean(PSD_PkRatio);
    averageFFTL = F(17);
    averagePSDL = F(18);
    FFTPeaks1 = F(19:22);
    PSDPeaks1 = F(23:26);
    b = F(27:30);
    if sum(b(1:2)) == 2 %FFTs match
        for i = 1:4
            if isequal(wLFFT, [i i i i])
                Y(2) = CLASS(i);
            end
        end
    end
    if sum(b(3:4))==2
        for i = 1:4
            if isequal(wLPSD, [i i i i])
                Y(3) = CLASS(i);
            end
        end
    else
        Y(3) = 0;
    end
    if sum(b) == 4
        Y(1) = -1;
        Y(4) = (Y(2) == Y(3)) && (Y(2)~=0); % IF Frequencies match up entirely. 
    end
%     FFTMatch
    unqwLFFT = unique(wLFFT);
    countFFT = histc(wLFFT, unqwLFFT);
    unqwPSD = unique(wLPSD);
    countPSD = histc(wLPSD, unqwPSD);
    if isequal(unqwLFFT, unqwPSD) && Y(4)
        Y(5) = 1;
    end
end

end %TreeClassifier Function

