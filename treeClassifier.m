function [ Y ] = treeClassifier( F, F2 )
%treeClassifier - For SSVEP Classification. 

CLASS = [10 12 15 16];

numFeatures = length(F2); %if L = 72
%% F1 %# OF FEATURES IS CONSTANT:
Y = zeros(7,1);
%short classifier
% unpack variables:
wLFFT = F(1:4); wLPSD = F(5:8);
FFT_PkRatio = F(9:12); PSD_PkRatio = F(13:16);
averagePkRatioFFT = mean(FFT_PkRatio);
averagePkRatioPSD = mean(PSD_PkRatio);
% averageFFTL = F(17); %May not be very useful. 
% averagePSDL = F(18);
FFTPeaks1 = F(19:22); % Locations of major peaks
PSDPeaks1 = F(23:26);
b = F(27:30);
A1 = 0;
A2 = 0;
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
%% F2
%Unpack shared data:
fft_sel_loc = zeros(4,4);
fft_sel_pks = zeros(4,4);
psd_sel_loc = zeros(4,4);
psd_sel_pks = zeros(4,4);
stft_sel_pks = zeros(1,4);
stft_sel_loc = zeros(1,4);

%COMPARE PEAKS
if numFeatures ==64
    for i = 1:4 %LOC BLOCKS:
        fft_sel_loc(i,:) = F2((i-1)*16+1:(i-1)*16+4);
        psd_sel_loc(i,:) = F2((i-1)*16+5:(i-1)*16+8);
        fft_sel_pks(i,:) = F2((i-1)*16+9:(i-1)*16+12);
        psd_sel_pks(i,:) = F2((i-1)*16+13:(i-1)*16+16);
    end
else %numFeatures = 72
    for i = 1:4 %LOC BLOCKS:
        fft_sel_loc(i,:) = F2((i-1)*18+1:(i-1)*18+4);
        psd_sel_loc(i,:) = F2((i-1)*18+5:(i-1)*18+8);
        fft_sel_pks(i,:) = F2((i-1)*18+9:(i-1)*18+12);
        psd_sel_pks(i,:) = F2((i-1)*18+13:(i-1)*18+16);
        stft_sel_loc(i) = F2((i-1)*18+17);
        stft_sel_pks(i) = F2((i-1)*18+18);
    end
    
end
%% Analysis F2
B1 = fft_sel_loc~=0;
sumB1 = sum(B1,2)
B2 = psd_sel_loc~=0;
sumB2 = sum(B2,2)
B1_1 = (sumB1 == 4);
B2_1 = (sumB2 == 4);
B_Short = B1_1 & B2_1; %B1_1 & B2_1; %IF ROWS IN FFT AND PSD ARE COMPLETE
if sum(B_Short) > 1
    %COMPARE PEAKS (SUM)
%     PkValuesToCompare = sum(fft_sel_pks(B_Short,:),2);
    PkVFFT = sum(fft_sel_pks,2);
    PkVPSD = sum(psd_sel_pks,2); 
    [~,I1] = max(PkVFFT);
    [~,I2] = max(PkVPSD);
    if(I1 == I2)
        A1 = CLASS(I1);
    else
        A1 = 0;
    end
elseif sum(B_Short) == 1
    A1 = CLASS(B_Short);
else
    A1 = 0;
end

if sum(A1)==1
    A1 = CLASS(A1);
end

if ~isempty(A1)
    if (A1(1) == Y(2)) && (Y(4) == 1)
        Y(1) = 1;
    end
else
    Y(1) = 0;
end

if numFeatures == 72
    B_L = (B_Short & ((stft_sel_loc')~=0));
    if sum(B_L)==1
        A2 = CLASS(B_L);
    else
        A2 = 0;
    end
    if isempty(A1) 
        A1 = 0;
    end
    if isempty(A2)
        A2 = 0;
    end
%     if (A1 == A2) && (Y(4)==1) && (A1 == Y(2))
    if (A1(1) == A2(1)) && (A1(1) ~= 0)
        Y(1) = 1;
    else
        Y(1) = 0;
    end
end
Y(6) = A1;
Y(7) = A2; 
% Y = [Y];

end %TreeClassifier Function

