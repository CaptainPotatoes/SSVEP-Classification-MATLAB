function [ Y ] = eogcfilt_b( X )
% Butterworth Order 3 [0.15-9.5Hz] for EOG.
X = X(:);
b = [0.00129953028712882,0,-0.00389859086138647,0,0.00389859086138647,0,-0.00129953028712900];
a = [1,-5.52855503341017,12.7497009825319,-15.7020189237831,10.8934835499658,-4.03693878488703,0.624328210166808];
Y = filtfilt(b,a,X);

end
