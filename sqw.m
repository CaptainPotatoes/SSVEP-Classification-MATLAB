%Sq. Wav.
clear;clc;close all;
T_period = [33 30 27 25];
f_SSVEP = 1000 * 1./(T_period.*2);
t = 0:0.001:1;
str = {'(LED 1) 15.15Hz', ...
    '(LED 2) 16.67Hz', ...
    '(LED 3) 18.52Hz', ...
    '(LED 4) 20.0Hz'};
for i=1:4
    figure(i); hold on;
    x = square(t*(2*pi*f_SSVEP(i)));%-4*i + 
    plot(t,x), ylim([-2 2])
    legend(str{i});
    
    figure(7); hold on;
    Y = fft(x);
    L = length(x);
    f = 1000*(0:(L/2))/L;
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(f,P1);
end
figure(1); %legend(str);%axis off
set(gca,'ytick',[]);
xlabel('Time (s)');
ylabel('15 Hz Stimulus Timing');
