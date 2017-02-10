function [ F ] = featureExtraction( classifyingWindow )
%featureExtraction Summary of this function goes here if I ever feel like
%writing one up.
T_mean = Wmean(classifyingWindow);
T_stdv = Wstd(classifyingWindow);
% T_pca  = Wpca1(classifyingWindow);
T_max = Wmax(classifyingWindow);
T_min = Wmin(classifyingWindow);
T_countmin_1 = sum(classifyingWindow<-1*10^-5 & classifyingWindow>(-1*10^-4),2);
T_countmin_2 = WCountMin(classifyingWindow,-1*10^-4);
T_countmax = WCountMax(classifyingWindow,8.5*10^-5);
% T_Integrate = WIntegrate(classifyingWindow);
T_Integrate = trapz(classifyingWindow);
NUMBER_FEATURES = 8;
% F = horzcat(T_mean, T_stdv, T_max, T_min, T_countmin_1, T_countmin_2, T_countmax, T_Integrate);
% F = [T_mean, T_stdv, T_pca, T_max, T_min, T_countmin, T_countmax, T_Integrate]; 
% FINE KNN:
% USE: 
F = horzcat( T_max, T_min, T_countmin_1, T_countmin_2, T_countmax);
end

