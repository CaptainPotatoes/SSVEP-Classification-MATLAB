function [ F ] = featureExtractionEOG2( samplesX )
%featureExtraction Summary of this function goes here if I ever feel like
%writing one up.
THRESHOLD1 = 2.85E-4;
THRESHOLD2 = 3.6E-3;
T_average_peak_distance = 0;
% figure(20);hold on; plot(samplesX);
% waveletCoef0 = cwt(samplesX,1:20,'haar');
% waveletCoef = cwt_haar(samplesX);
% Scalogram = abs(waveletCoef.*waveletCoef);
% E = sum(Scalogram(:)); %Energy Level.
IntegralEst = trapz(samplesX);
[Max,Imax] = findpeaks(diff(samplesX), 'SORTSTR','descend','NPEAKS',1);
[Min,Imin] = findpeaks(-diff(samplesX),'SORTSTR','descend','NPEAKS',1);
Amplitude = abs(Max-Min);
Velocity = Amplitude/(Imax - Imin);
Mean = (1/size(samplesX,2))*sum(samplesX);
T_stdv = Wstd(samplesX);
T_max = Wmax(samplesX);
T_min = Wmin(samplesX);
T_Integrate = trapz(samplesX);
peaks = [];
T_findpeaks_distX=[];
[peaks, loc] = findpeaks(samplesX, 'MinPeakHeight', 1E-3);
if isempty(peaks)
    T_count_findpeaks = 0;
    T_findpeaks_distX = 0;
else
    T_count_findpeaks = length(peaks);
    lp = zeros(length(peaks),1);
    if length(peaks)>1
        for i = 2:length(peaks)
            lp(i-1) = loc(i)-loc(i-1);
        end
        T_average_peak_distance = mean(lp); 
        T_findpeaks_distX = loc(end) - loc(1); %TODO: TAKE AVG, NOT MAX-MIN
    else
        T_findpeaks_distX = 0;
    end
end
%%? Threshold Stuff:?
T_countmin_1 = sum(samplesX<THRESHOLD2 & samplesX>(THRESHOLD1),2);
T_countmin_2 = WCountMin(samplesX,THRESHOLD1);
% T_countmax = WCountMax(samplesX,THRESHOLD2);

% F = horzcat(E,IntegralEst, Velocity, Amplitude, Mean, T_stdv, T_max, T_min, T_Integrate, T_count_findpeaks, T_findpeaks_distX, T_average_peak_distance, T_countmax, T_countmin_1, T_countmin_2);
F = horzcat(Amplitude, Velocity, Mean, T_stdv, T_max, T_min, T_Integrate, T_count_findpeaks, T_findpeaks_distX, T_average_peak_distance, T_countmin_1, T_countmin_2);
% hold off;
end

