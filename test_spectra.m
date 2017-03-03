clear all; clc; close all
% model parameters for test
r = 0.99; % Radius of the first pole
kk = 1:5;
rr = r.^kk; % all poles' radius
thetha = pi/7; % First pole angle
thethar = thetha.*kk; % all poles' angle
rN = 14;
rM = 4;
Burnin = 10000;
x = randn(1,Burnin+2^(rN+1)); % white noise
% para gerar espectro te?rico do sinal de teste
xa = 1;
for i=1:5
% 'xa' will have the coefficients for the theoretical spectra
xa = conv([1 -2*rr(i)*cos(thethar(i)) rr(i)^2],xa);
x = filter(1,[1 -2*rr(i)*cos(thethar(i)) rr(i)^2],x); % filter the signal
% to then generate 5 tonal peaks
end
% burn the transient
x = x(Burnin+1:end);
%%
linewidth = 1.3;
x = x(:);
disp(['signal variance: ' num2str(var(x))])
figure, hold on;
[ft,w]=freqz(1,xa,1024);
plot(w*.5/pi,20*log10(abs(ft)),'r.','LineWidth',linewidth)
wind = [1024 512 256 128];
for i=1:length(wind)
[spec,f]=welch_estimator(x,1,hann(wind(i)));
plot(f,10*log10(abs(reshape(spec/2,length(f),1))),'LineWidth',linewidth)
disp(['PSD integral w/ L=' num2str(wind(i)) ': ' num2str(sum(spec)*f(2))])
end
legend( 'theoretical spec',...
['Window ' num2str(wind(1))],...
['Window ' num2str(wind(2))],...
['Window ' num2str(wind(3))],...
['Window ' num2str(wind(4))])
xlabel('Normalized frequency'); ylabel('PSD [dB]');
title('Power spectra test');
grid
set(gca, 'FontSize', 13)
