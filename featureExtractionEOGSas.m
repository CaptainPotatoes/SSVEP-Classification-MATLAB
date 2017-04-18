function F = eogFeaturesSas(ClassifyingWindow1,ClassifyingWindow2,ClassifyingWindow3,ClassifyingWindow4,ClassifyingWindow5,ClassifyingWindow6)
%Here is two inputs instead of;
%CW1 is the filtered signal 
%CW2 is the difference of the filtered signal (filter the filtered and creatge a two point differentiation)

waveletcoef = cwt(ClassifyingWindow1,1:20,'haar');% Obtains coefficients for wavelet transform
S = abs(waveletcoef.*waveletcoef); %Obtain Scalogram
E = sum(S(:));%Obtain Energy Level
TRAPZ = trapz(ClassifyingWindow1);
[ Maxima, MaxIdx] = findpeaks(diff(ClassifyingWindow1),'SORTSTR','descend','NPEAKS',1);
[ Minima, MinIdx] = findpeaks(-diff(ClassifyingWindow1),'SORTSTR','descend','NPEAKS',1);
Minima = - Minima;
Amplitude = abs( Maxima -  Minima);
Velocity =  Amplitude/( MaxIdx- MinIdx);
Mean = (1/size(ClassifyingWindow1,2))*sum(ClassifyingWindow1);
CUMTRAPZ = sum(cumtrapz(ClassifyingWindow1));
COV = cov(ClassifyingWindow1);
RMS = rms(ClassifyingWindow1);
PEAKTOPEAK = peak2peak(ClassifyingWindow1);
PEAKTORMS = peak2rms(ClassifyingWindow1);
% MAX = max(ClassifyingWindow1);
% MIN = min(ClassifyingWindow1);



waveletcoef2 = cwt(ClassifyingWindow2,1:20,'haar');% Obtains coefficients for wavelet transform
S2 = abs(waveletcoef2.*waveletcoef2); %Obtain Scalogram
E2 = sum(S2(:));%Obtain Energy Level
TRAPZ2 = trapz(ClassifyingWindow2);
[ Maxima2, MaxIdx2] = findpeaks(diff(ClassifyingWindow2),'SORTSTR','descend','NPEAKS',1);
[ Minima2, MinIdx2] = findpeaks(-diff(ClassifyingWindow2),'SORTSTR','descend','NPEAKS',1);
Minima2 = - Minima2;
Amplitude2 = abs( Maxima2 -  Minima2);
Velocity2 =  Amplitude2/( MaxIdx2- MinIdx2);
Mean2 = (1/size(ClassifyingWindow2,2))*sum(ClassifyingWindow2);
CUMTRAPZ2 = sum(cumtrapz(ClassifyingWindow2));
COV2 = cov(ClassifyingWindow2);
RMS2 = rms(ClassifyingWindow2);
PEAKTOPEAK2 = peak2peak(ClassifyingWindow2);
PEAKTORMS2 = peak2rms(ClassifyingWindow2);

waveletcoef3 = cwt(ClassifyingWindow3,1:20,'haar');% Obtains coefficients for wavelet transform
S3 = abs(waveletcoef3.*waveletcoef3); %Obtain Scalogram
E3 = sum(S3(:));%Obtain Energy Level
TRAPZ3 = trapz(ClassifyingWindow3);
[ Maxima3, MaxIdx3] = findpeaks(diff(ClassifyingWindow3),'SORTSTR','descend','NPEAKS',1);
[ Minima3, MinIdx3] = findpeaks(-diff(ClassifyingWindow3),'SORTSTR','descend','NPEAKS',1);
Minima3 = - Minima3;
Amplitude3 = abs( Maxima3 -  Minima3);
Velocity3 =  Amplitude3/( MaxIdx3- MinIdx3);
Mean3 = (1/size(ClassifyingWindow3,2))*sum(ClassifyingWindow3);
CUMTRAPZ3 = sum(cumtrapz(ClassifyingWindow3));
COV3 = cov(ClassifyingWindow3);
RMS3 = rms(ClassifyingWindow3);
PEAKTOPEAK3 = peak2peak(ClassifyingWindow3);
PEAKTORMS3 = peak2rms(ClassifyingWindow3);
% MAX = max(ClassifyingWindow1);
% MIN = min(ClassifyingWindow1);



waveletcoef4 = cwt(ClassifyingWindow4,1:20,'haar');% Obtains coefficients for wavelet transform
S4 = abs(waveletcoef4.*waveletcoef4); %Obtain Scalogram
E4 = sum(S4(:));%Obtain Energy Level
TRAPZ4 = trapz(ClassifyingWindow4);
[ Maxima4, MaxIdx4] = findpeaks(diff(ClassifyingWindow4),'SORTSTR','descend','NPEAKS',1);
[ Minima4, MinIdx4] = findpeaks(-diff(ClassifyingWindow4),'SORTSTR','descend','NPEAKS',1);
Minima4 = - Minima4;
Amplitude4 = abs( Maxima4 -  Minima4);
Velocity4 =  Amplitude4/( MaxIdx4- MinIdx4);
Mean4 = (1/size(ClassifyingWindow4,2))*sum(ClassifyingWindow4);
CUMTRAPZ4 = sum(cumtrapz(ClassifyingWindow4));
COV4 = cov(ClassifyingWindow4);
RMS4 = rms(ClassifyingWindow4);
PEAKTOPEAK4 = peak2peak(ClassifyingWindow4);
PEAKTORMS4 = peak2rms(ClassifyingWindow4);


waveletcoef5 = cwt(ClassifyingWindow5,1:20,'haar');% Obtains coefficients for wavelet transform
S5 = abs(waveletcoef5.*waveletcoef5); %Obtain Scalogram
E5 = sum(S5(:));%Obtain Energy Level
TRAPZ5 = trapz(ClassifyingWindow5);
[ Maxima5, MaxIdx5] = findpeaks(diff(ClassifyingWindow5),'SORTSTR','descend','NPEAKS',1);
[ Minima5, MinIdx5] = findpeaks(-diff(ClassifyingWindow5),'SORTSTR','descend','NPEAKS',1);
Minima5 = - Minima5;
Amplitude5 = abs( Maxima5 -  Minima5);
Velocity5 =  Amplitude5/( MaxIdx5- MinIdx5);
Mean5 = (1/size(ClassifyingWindow5,2))*sum(ClassifyingWindow5);
CUMTRAPZ5 = sum(cumtrapz(ClassifyingWindow5));
COV5 = cov(ClassifyingWindow5);
RMS5 = rms(ClassifyingWindow5);
PEAKTOPEAK5 = peak2peak(ClassifyingWindow5);
PEAKTORMS5 = peak2rms(ClassifyingWindow5);
% MAX = max(ClassifyingWindow1);
% MIN = min(ClassifyingWindow1);



waveletcoef6 = cwt(ClassifyingWindow6,1:20,'haar');% Obtains coefficients for wavelet transform
S6 = abs(waveletcoef6.*waveletcoef6); %Obtain Scalogram
E6 = sum(S6(:));%Obtain Energy Level
TRAPZ6 = trapz(ClassifyingWindow6);
[ Maxima6, MaxIdx6] = findpeaks(diff(ClassifyingWindow6),'SORTSTR','descend','NPEAKS',1);
[ Minima6, MinIdx6] = findpeaks(-diff(ClassifyingWindow6),'SORTSTR','descend','NPEAKS',1);
Minima6 = - Minima6;
Amplitude6 = abs( Maxima6 -  Minima6);
Velocity6 =  Amplitude6/( MaxIdx6- MinIdx6);
Mean6 = (1/size(ClassifyingWindow6,2))*sum(ClassifyingWindow6);
CUMTRAPZ6 = sum(cumtrapz(ClassifyingWindow6));
COV6 = cov(ClassifyingWindow6);
RMS6 = rms(ClassifyingWindow6);
PEAKTOPEAK6 = peak2peak(ClassifyingWindow6);
PEAKTORMS6 = peak2rms(ClassifyingWindow6);
% MAX2 = max(ClassifyingWindow2);
% MIN2 = min(ClassifyingWindow2);
%             Feature = [ TRAPZ, Velocity,  Amplitude,  E, Mean];
F = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ,COV,RMS,PEAKTOPEAK,PEAKTORMS,TRAPZ2, Velocity2 Amplitude2, E2, Mean2, CUMTRAPZ2,COV2,RMS2,PEAKTOPEAK2,PEAKTORMS2 ...
     TRAPZ3,Velocity3,Amplitude3,E3,Mean3,CUMTRAPZ3,COV3,RMS3,PEAKTOPEAK3,PEAKTORMS3,TRAPZ4, Velocity4 Amplitude4, E4, Mean4, CUMTRAPZ4,COV4,RMS4,PEAKTOPEAK4,PEAKTORMS4...
     TRAPZ5,Velocity5,Amplitude5,E5,Mean5,CUMTRAPZ5,COV5,RMS5,PEAKTOPEAK5,PEAKTORMS5,TRAPZ6, Velocity6 Amplitude6, E6, Mean6, CUMTRAPZ6,COV6,RMS6,PEAKTOPEAK6,PEAKTORMS6];
% Feature = [ CUMTRAPZ,COV,Velocity,Amplitude,Mean,E];
% Feature = [ CUMTRAPZ];
% if NumberofFeatures == 1
% Feature = [ TRAPZ,TRAPZ2];
% elseif NumberofFeatures == 2
% Feature = [ TRAPZ,Velocity,TRAPZ2, Velocity2]; 
% elseif NumberofFeatures == 3
% Feature = [ TRAPZ,Velocity,Amplitude,TRAPZ2, Velocity2, Amplitude2];   
% elseif NumberofFeatures == 4
% Feature = [ TRAPZ,Velocity,Amplitude,E,TRAPZ2, Velocity2 Amplitude2, E2];   
% elseif NumberofFeatures == 5
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,TRAPZ2, Velocity2, Amplitude2, E2, Mean2];  
% elseif NumberofFeatures == 6
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ,TRAPZ2, Velocity2 Amplitude2, E2, Mean2, CUMTRAPZ2];
% elseif NumberofFeatures == 7
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ,COV, Velocity2, Amplitude2, E2, Mean2, CUMTRAPZ2,COV2];
% elseif NumberofFeatures == 8
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ,COV,RMS,TRAPZ2, Velocity2 Amplitude2, E2, Mean2, CUMTRAPZ2,COV2,RMS2];    
% elseif NumberofFeatures == 9    
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ,COV,RMS,PEAKTOPEAK,TRAPZ2, Velocity2, Amplitude2, E2, Mean2, CUMTRAPZ2,COV2,RMS2,PEAKTOPEAK2];
% elseif NumberofFeatures == 10 
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ,COV,RMS,PEAKTOPEAK,PEAKTORMS,TRAPZ2, Velocity2, Amplitude2, E2, Mean2, CUMTRAPZ2,COV2,RMS2,PEAKTOPEAK2,PEAKTORMS2];   
% end

