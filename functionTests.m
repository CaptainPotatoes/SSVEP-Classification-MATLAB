clear;clc;close all;
load('mssvep_16.6_3.mat');

ch1 = Trial{1}(1:end,1); %ignore last second
ch1_f = eeg_h_custom(ch1, 250, [8 20], 5);
ch1_f_n = scaleAbs(ch1_f);
figure(1);
% hold on;
plot(ch1_f);
figure(2);
plot(ch1_f_n);
% hold off;

