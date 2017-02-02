%% Import
clear all;clc;close all;
% Open File:
channelNames = {'Fp1' 'Fp2'};
[filename,pathname,~]  = uigetfile('*.csv','Select Data to Upload');
recording1 = csvread([pathname filename]);

Fs = 250;
% fp1_recording1 = recording1(:,1);
%%TEMP 
fp1_recording1 = recording1(:,1);
% fp1_recording1 = BPSignal{1};
fp2_recording1 = recording1(:,2);
% fp2_recording1 = BPSignal{2};

h=1/Fs;
L = size(fp1_recording1,1);
t = 0:h:L/Fs-h;
figure(10);
plot(fp1_recording1(1:500));
filtEEG = eeg_bandpass(fp1_recording1(1:500), Fs);
figure(11);
plot(filtEEG);
figure(1);
subplot(2,1,1);
plot(t,fp1_recording1);
% hold on;
subplot(2,1,2);
plot(t,fp2_recording1);
% hold off;

start = 1;
seconds = 5;
length = seconds*Fs+1;   
% window_seconds = length/Fs
n_cells = floor(size(fp1_recording1,1)/length)
for i=0:n_cells-1
    fp1{i+1}=fp1_recording1(start+i*length:start+length*(i+1));
end
for i=0:n_cells-1
    fp2{i+1}=fp2_recording1(start+i*length:start+length*(i+1));
end
%% Spectrogram (Full Recording)
    % Fp1:
%     tic
    %Filter Fp1:
    fp1_r_filtered = eeg_h_fcn(fp1_recording1(1:1250),Fs);
%     fp1_r_filtered = eeg_bandpass(fp1_recording1,Fs);
%     fp2_r_filtered = eeg_h_fcn(fp2_recording1,Fs);
    n=1;
    %Formof: [s,f,t,p] = spectrogram(x,window,noverlap,f,fs) 
        % f = cyclical frequencies
    startFromHz = 0;
    upToHz = 240;
%     [S{n},Fspect{n},T{n},P{n}] = spectrogram(fp1_r_filtered(:,n),5*Fs,4*Fs,10*Fs,Fs);
    [S{n},Fspect{n},T{n},P{n}] = spectrogram(fp1_r_filtered, Fs,0.5*Fs,10*Fs,Fs);
    figure
    imagesc( T{n}, ...
             Fspect{n}(Fspect{n}<upToHz & Fspect{n}>startFromHz), ...
             10*log10(P{n}(Fspect{n}<upToHz & Fspect{n}>startFromHz,:)) )
%     toc
    set(gca,'YDir','normal')
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    cb = colorbar;
    ylabel(cb, 'Power (db)')
    colormap(jet)
    title(['Channel ' channelNames{n}], 'FontSize', 14)
    %Fp2:
%     n=2;
%     %Formof: [s,f,t,p] = spectrogram(x,window,noverlap,f,fs) 
%         % f = cyclical frequencies
%     upToHz = 50;
%     [S{n},Fspect{n},T{n},P{n}] = spectrogram(fp2_r_filtered(:,1),5*Fs,4*Fs,10*Fs,Fs);
%     figure
%     imagesc( T{n}, ...
%              Fspect{n}(Fspect{n}<upToHz), ...
%              10*log10(P{n}(Fspect{n}<upToHz,:)) )
% %     toc
%     set(gca,'YDir','normal')
%     ylabel('Frequency (Hz)')
%     xlabel('Time (s)')
%     cb = colorbar;
%     ylabel(cb, 'Power (db)')
%     colormap(jet)
%     title(['Channel ' channelNames{n}], 'FontSize', 14)
% % clear recording1;
%% Averaging Filter:
% window = 100;n=1;
% for i = 1:length(fp1_recording1)
%     if i >= window && i < length(fp1_recording1)-window
%         meandata(n,i) = mean(fp1_recording1(i-window/2:i+window/2,n));
%     elseif i < window
%         meandata(n,i) = fp1_recording1(i,n);
%     else
%         meandata(n,i) = mean(fp1_recording1...
%             (i-window/2:length(fp1_recording1(:,2))));
%     end
% end
% avgdata = meandata';
% shiftedData(:,n) = fp1_recording1(:,n) - avgdata(:,n);
% figure
% % subplot (2,1,1);
% plot(shiftedData(:,1));

%% Filter (BW, O3, )
L2 = size(fp1{1},1);
for i=1:n_cells
%     figure
    %%fp1_f{i} = eegFir(fp1{i}, Fs);
    fp1_f{i} = eeg_h_fcn(fp1{i}, Fs);
    fp2_f{i} = eeg_h_fcn(fp2{i}, Fs);
    t2 = 0:h:L2/Fs-h;
%     subplot(2,1,1);
%     plot(t2,fp1_f{i});
%     subplot(2,1,2);
%     plot(t2,fp2_f{i});
end

%% Spectrogram of Shorter Signals:
seconds
    tic
    i=5;
    n=1;
    dBmax = 100;
    [S{n},Fspect{n},T{n},P{n}] = spectrogram(fp1_f{i},Fs,0.5*Fs,2*Fs,Fs);
    figure
    
    imagesc(T{n},Fspect{n}(Fspect{n}<dBmax),10*log10(P{n}(Fspect{n}<dBmax,:)))
    toc
    set(gca,'YDir','normal')
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    cb = colorbar;
    ylabel(cb, 'Power (db)')
    colormap(jet)
    title(['Channel ' channelNames{n}], 'FontSize', 14)


%% FFT (maybe implement a smoothing filter and see how it affects measurements)
% figure(3)
% hold on;
for i=1:n_cells
    figure
    fp1_fft{i} = fft(fp1_f{i});
    P2 = abs(fp1_fft{i}/length);
    P1 = P2(1:length/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f=Fs*(0:(length/2))/length;
    plot(f,P1),xlim([0,100]),ylim([0, max(P1)]);
    title('FFT of SSVEP Data [Fp1]');
    xlabel('f (Hz)');
    ylabel('Magnitude (\muV)');
    %maxP1_40 = max(P1(116:127))
%     peakP1_46(i) = max(P1(133:145));
%     test1{i} = maxP1>0.060;
end

% minMax = [min(peakP1_46) max(peakP1_46)]
% rangeMaxP1_46 = range(peakP1_46)
% for i=1:n_cells
%     fprintf('i = ');
%     if test1{i}
%         fprintf('40Hz\r\n');
%     else
%         fprintf('Not 40Hz\r\n');
%     end
% end
% hold off;

% Normal FFT:
% figure(4);
% hold on;
% fft_n=2^nextpow2(size(fp1_f{1},1));
% for i=1:n_cells
%     fp1_fft2{i} = fft(fp1{i}, fft_n); % using next power of two
%     P = abs(fp1_fft2{i}/fft_n);
%     f_2=Fs*(0:(fft_n/2))/fft_n;
%     plot(f_2,P(1:fft_n/2+1)),xlim([1,70]);
% end
% hold off;
% peakFFT = max(P1);
% P1_2 = P1(5:size(P1,1));
% peakFFT2 = max(P1(5:size(P1,1)));
%% Display peak frequency:
% Need to make this more complex. Need to find multiple peak values
    % Looking for 35.7Hz
% [c1, c2] = find(ismember(P1, max(P1(:))))
% For 2s:
m = 2;
[b1, b2] = find(ismember(P1, max(P1(seconds*35:45*seconds))))
[c1, c2] = find(ismember(P1, max(P1(seconds*68:seconds*75))))
[d1, d2] = find(ismember(P1, max(P1(seconds*76:seconds*86))))
[e1, e2] = find(ismember(P1, max(P1(seconds*87:seconds*99))))
% FOR 1s:
% [c1, c2] = find(ismember(P1, max(P1(69:75))))
% [d1, d2] = find(ismember(P1, max(P1(76:86)))).
peak_frequency = f(1, b1)
size = P1(b1)
peak_frequency = f(1, c1)
size = P1(c1)
peak_frequency2 = f(1, d1)
size2 = P1(d1)
peak_frequency3 = f(1, e1)
size3 = P1(e1)
peakdiff = size3/size2
if abs(peak_frequency-71.43)<4%Hz
    fprintf('Success \n');
else
    fprintf('Fail \n');
end
%%
clear all;clc;close all;
% load('Trial2Data.mat')
tD = cell(1);
% load('BaseLineData1')
load('BaseLineData2')
fp1d = Trial{1}(1:end-250,1); %ignore last second
fp2d = Trial{2}(1:end-250,1);
h=1/250;
t=0:h:length(fp1d)/250-h;
markers = tD{1};
% figure
% hold on;
%     plot(fp1d);
%     plot(fp2d);
%     for i=1:length(markers)
%         text(markers(i,1), fp1d(markers(i,1)), num2str(markers(i,2)));
%     end
% hold off;
% Filter Entirety
fp1dfilt = eog_h_fcn(fp1d,250);
fp2dfilt = eog_h_fcn(fp2d,250);
figure
plot(t,fp1dfilt,'color','g'),ylim([-8e-4,8e-4]);
for i=1:length(markers)
    text(t(markers(i,1)), fp1dfilt(markers(i,1)), num2str(markers(i,2)));
end

figure
plot(t,fp2dfilt,'color','r'),ylim([-8e-4,8e-4]);
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
Features(:,numFeatures+1) = [Class];

%% Sas Method
%separate into wndows then filter & extract features
winLen = 250;
Window = cell(floor(length(fp1d)/winLen),1);
% Features = cell(floor(length(fp1d)/winLen),1);
for i = 1:floor(length(fp1d)/winLen)
    start = (i-1)*winLen+1;
%     Window{i} = fp1d(start:start+winLen-1);
    Window{i} = eog_h_fcn(fp1d(start:start+winLen-1),250);
%     Features{i,2} = Feature(Window{i});
    Features(i,:) = [Feature(Window{i}),-1];
%     CLASS:
    for j=1:length(tD{1})
        %if window contains number, assign class, else zero
        if sum((start:start+winLen-1)==tD{1}(j,1))
%             Features{i,1} = tD{1}(j,2);
            Features(i,11) = tD{1}(j,2);
            break;
        else
%             Features{i,1} = 0;
            Features(i,11) = 0;
        end
    end
end