%% Load, plot and annotate.
clear all;clc;close all;
% load('Trial2Data.mat')
% tD = cell(1);
% load('BaseLineData1')
% load('BaseLineData2')
% load('Trial_DB1');
load('Trial_Marc_2.mat');
fp1d = Trial{1}(1:end-250,1); %ignore last second
fp2d = Trial{2}(1:end-250,1);
h=1/250;
t=0:h:length(fp1d)/250-h;
markers = tD{1};

% Filter Entirety
fp1dfilt = eog_h_fcn(fp1d,250);
fp2dfilt = eog_h_fcn(fp2d,250);
figure
plot(t,fp1dfilt,'color','g')%,ylim([-8e-4,8e-4]);
for i=1:length(markers)
    text(t(markers(i,1)), fp1dfilt(markers(i,1)), num2str(markers(i,2)));
end

figure
plot(t,fp2dfilt,'color','r')%,ylim([-8e-4,8e-4]);
for i=1:length(markers)
    text(t(markers(i,1)), fp2dfilt(markers(i,1)), num2str(markers(i,2)));
end
%% Gandhi Method
winLen = 500;
Window = cell(floor(length(fp1d)/winLen),1);
% numFeatures = 5;
% Features = zeros(floor(length(fp1d)/winLen),numFeatures+1);
% Features = zeros(120,9);
% Features = cell(floor(length(fp1d)/winLen),1);
hold all;
for i = 1:floor(length(fp1d)/winLen)
    start = (i-1)*winLen+1;
    Window{i} = eog_h_fcn(fp1d(start:start+winLen-1),250);
%     xalloc = ((i-1)*winLen+1):(i*winLen);
%     yalloc = 1:numFeatures;
%     Features(xalloc,yalloc) = featureExtraction(Window{i});
    [numFeatures Features(i,:)] = featureExtraction(Window{i}');
%     CLASS:
%     Class = zeros (xalloc, 1);
%     for j=((i-1)*winLen+1):(i*winLen)
%         if sum(j>tD{1}(:,1)-250 & j<tD{1}(:,1)+250)
%             logMat = j>tD{1}(:,1)-250 & j<tD{1}(:,1)+250;
%             findV = find(logMat,1,'first');
%             Class(j,1) = tD{1}(findV,2);
%         else
%             Class(j,1) = 0;
%         end
%     end
    for j=1:length(tD{1})
        %if window contains number, assign class, else zero
        if sum((start:start+winLen-1)==tD{1}(j,1))
%             Features{i,1} = tD{1}(j,2);
            Class(i,1) = tD{1}(j,2);
            break;
        else
%             Features{i,1} = 0;
            Class(i,1) = 0;
        end
    end
end
Features(:,numFeatures+1) = Class;
%% Shifting Window Method
clear;clc;close all;
% load('Trial_DB1');
load('Trial_Marc_3.mat')
fp1 = Trial{1}(1:end-250,1); %ignore last second
fp2 = Trial{2}(1:end-250,1);
h=1/250;
t=0:h:length(fp1)/250-h;
% markers = tD{1};
markers = trainingData{1};
Fs = SamplingRate; 
% Fs = 250;
fp1filt = eog_h_fcn(fp1,250);
fp2filt = eog_h_fcn(fp2,250);
% figure(1)
% hold on;
% plot(t,fp1filt,'color','r')%,ylim([-8e-4,8e-4]);
% 
% plot(t,fp2filt,'color','y')%,ylim([-8e-4,8e-4]);
% for i=1:length(markers)
%     text(t(markers(i,1)), fp2filt(markers(i,1)), num2str(markers(i,2)));
% end
% hold off

figure(2)
hold on;
plot(fp1filt,'color','r')%,ylim([-8e-4,8e-4]);

plot(fp2filt,'color','y')%,ylim([-8e-4,8e-4]);
for i=1:length(markers)
    text(markers(i,1), fp2filt(markers(i,1)), num2str(markers(i,2)));
end
hold off
%% Set classes:
% close all;
seconds = 2; %2 second window
winLen = seconds*Fs; 
winFraction = 2;%2.5; %1/5 of a second
winShift = floor(Fs/winFraction); 
dataLimit = floor((length(fp1)-winLen)/winLen);
start = 1;
Window = cell( seconds*winFraction*dataLimit - 1, 2);
% Window = cell( 10*dataLimit - 1, 2);
%Todo: preallocate F_fp1/2
assignedClass = zeros( seconds*winFraction*dataLimit - 1, 1);

figure(1); 
for i = 1 : seconds*winFraction*dataLimit
    start = 1 + winShift*(i-1);
    winEnd = start + winLen-1;
    fprintf('Current index = [%d to %d]\r\n',start, winEnd);
    Window{i,1} = fp1( start : start + winLen-1 );              % set values:
    Window{i,2} = fp2( start : start + winLen-1 );
    fp1f = eog_h_fcn( Window{i,1}, Fs); 
    fp2f = eog_h_fcn( Window{i,2}, Fs);
    hold on;
    plot(fp1f),ylim([-2.5E-4 2.5E-4]);    % Plot filtered Data. 
	plot(fp2f),ylim([-2.5E-4 2.5E-4]);
    hold off;
    assignedClass(i,1) = input('Enter an integer value!\n');
    F_fp1(i,:) = featureExtraction( fp1f' );
    F_fp2(i,:) = featureExtraction( fp2f' );
    clf(1);
end
%todo: horzcat [ Features ; assignedClass ] 
% ALLFEATURES = [F_fp1 F_fp2 assignedClass];

