function [ Y ] = eogcfilt( X )
%EOGCFILT EOG Filter for conversion to C. All inputs must be constant.
% Vectorize:
X = X(:);
% Sampling Frequency = 250;
%BW for 10Hz upper bound, Order of 3.
b = [0.00156701035058832,0.00470103105176495,0.00470103105176495,0.00156701035058832];
a = [1,-2.49860834469118,2.11525412700316,-0.604109699507275];
Z = filtfilt(b,a,X);
%BW filt for 2Hz lower bound, Order of 3:
d = [0.950971887923409,-2.85291566377023,2.85291566377023,-0.950971887923409];
c = [1,-2.89947959461186,2.80394797738300,-0.904347531392409];
Y = filtfilt(d,c,Z);
end

