clear;clc;
load('tXtY.mat');
load('meog_t1.mat');
fp1 = Trial{1}(1:end-250,1); %ignore last second
fp2 = Trial{2}(1:end-250,1);
fpz = Trial{3}(1:end-250,1);
eyeR = Trial{4}(1:end-250,1);
Fs = SamplingRate; 
seconds = 1; %2 second window
winLen = seconds*Fs; 
winFraction = 2;%2.5; %1/4 of a second
winShift = floor(Fs/winFraction); 
dataLimit = floor((length(fp1)-winLen)/winLen);
start = 1;
numCh = 4;
Window = cell( seconds*winFraction*dataLimit - 1, numCh);
assignedClass = zeros( seconds*winFraction*dataLimit - 1, 1);
% figNum = 2;
% fH = figure(figNum); 
% set(fH, 'Position', [100, 100, 1200, 900]);
cont = [];
for i = 1 : seconds*winFraction*dataLimit
    start = 1 + winShift*(i-1);
    winEnd = start + winLen-1;
    fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    Window{i,1} = fp1( start : start + winLen-1 );              % set values:
    Window{i,2} = fp2( start : start + winLen-1 );
    Window{i,3} = fpz( start : start + winLen-1 );
    Window{i,4} = eyeR( start : start + winLen-1 );
    fp1f = eog_h_fcn( Window{i,1}, Fs); 
    fp2f = eog_h_fcn( Window{i,2}, Fs);
    fpzf = eog_h_fcn( Window{i,3}, Fs);
    eyeRf = eog_h_fcn( Window{i,4}, Fs);
%     hold on;
%     plot(fp1f),ylim([-2.5E-4 2.5E-4]);    % Plot filtered Data. 
%     plot(fp2f),ylim([-2.5E-4 2.5E-4]);
%     plot(fpzf),ylim([-2.5E-4 2.5E-4]);
%     plot(eyeRf),ylim([-2.5E-4 2.5E-4]);
%     hold off;
    Y(i) = fullHybridClassifier(Window{i,1}, Window{i,2}, Window{i,3}, Window{i,4}, tX, tY);
%     if isempty(cont)
%         cont = input('Continue? \n');
%     end
%     clf(2);
end
%% **********
% start = 10;
% for i = 1:1
%     data1(:,i) = Window{start+i,1};
% 
%     data2(:,i) = Window{start+i,2};
% 
%     data3(:,i) = Window{start+i,3};
% 
%     data4(:,i) = Window{start+i,4};
% end
% Y = fullHybridClassifier(data1, data2, data3, data4, tX, tY)';