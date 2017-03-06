function [ F ] = featureExtractionEOG( classifyingWindow )
%featureExtraction Summary of this function goes here if I ever feel like
%writing one up.
% classifyingWindow = classifyingWindow(:);
T_mean = Wmean(classifyingWindow);
T_stdv = Wstd(classifyingWindow);
T_max = Wmax(classifyingWindow);
T_min = Wmin(classifyingWindow);
T_countmin_1 = sum(classifyingWindow<-1*10^-5 & classifyingWindow>(-1*10^-4),2);
T_countmin_2 = WCountMin(classifyingWindow,-1*10^-4);
T_countmax = WCountMax(classifyingWindow,8.5*10^-5);
T_Integrate = trapz(classifyingWindow);
peaks = [];
T_findpeaks_distX=[];
[peaks, loc] = findpeaks(classifyingWindow, 'MinPeakHeight', 0.65E-4);
if isempty(peaks)
    T_count_findpeaks = 0;
    T_findpeaks_distX = 0;
else
    T_count_findpeaks = length(peaks);
    if length(peaks)>1
        T_findpeaks_distX = loc(end) - loc(1); %TODO: TAKE AVG, NOT MAX-MIN
    else
        T_findpeaks_distX = 0;
    end
end
% F = horzcat(T_mean, T_stdv, T_max, T_min, T_countmin_1, T_countmin_2, T_countmax, T_Integrate);
F = horzcat(T_mean, T_stdv, T_max, T_min, T_countmin_1, T_countmin_2, T_countmax, T_Integrate, T_count_findpeaks, T_findpeaks_distX);
end

