function Feature = Feature(ClassifyingWindow1)


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

%             Feature = [ TRAPZ, Velocity,  Amplitude,  E, Mean];
Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ,COV,RMS,PEAKTOPEAK,PEAKTORMS];

% Feature = [ CUMTRAPZ,COV,Velocity,Amplitude,Mean,E];
% Feature = [ CUMTRAPZ];
% if NumberofFeatures == 1
% Feature = [TRAPZ];    
% elseif NumberofFeatures == 2
% Feature = [ TRAPZ,Velocity];   
% elseif NumberofFeatures == 3
% Feature = [ TRAPZ,Velocity,Amplitude];   
% elseif NumberofFeatures == 4
% Feature = [ TRAPZ,Velocity,Amplitude,E];    
% elseif NumberofFeatures == 5
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean];    
% elseif NumberofFeatures == 6
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ];    
% elseif NumberofFeatures == 7
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ,COV];
% elseif NumberofFeatures == 8
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ,COV,RMS];    
% elseif NumberofFeatures == 9    
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ,COV,RMS,PEAKTOPEAK];
% elseif NumberofFeatures == 10 
% Feature = [ TRAPZ,Velocity,Amplitude,E,Mean,CUMTRAPZ,COV,RMS,PEAKTOPEAK,PEAKTORMS];    
% end

