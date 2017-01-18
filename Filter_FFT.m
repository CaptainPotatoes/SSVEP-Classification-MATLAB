%% Import
clear all;clc;close all;
% Open File:
channelNames = {'Fp1' 'Fp2'};
[filename,pathname,~]  = uigetfile('*.csv','Select Data to Upload');
recording1 = csvread([pathname filename]);

Fs = 500;
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
    fp1_r_filtered = eeg_bandpass(fp1_recording1,Fs);
%     fp1_r_filtered = eeg_bandpass(fp1_recording1,Fs);
%     fp2_r_filtered = eeg_h_fcn(fp2_recording1,Fs);
    n=1;
    %Formof: [s,f,t,p] = spectrogram(x,window,noverlap,f,fs) 
        % f = cyclical frequencies
    startFromHz = 0;
    upToHz = 240;
%     [S{n},Fspect{n},T{n},P{n}] = spectrogram(fp1_r_filtered(:,n),5*Fs,4*Fs,10*Fs,Fs);
    [S{n},Fspect{n},T{n},P{n}] = spectrogram(fp1_recording1,5*Fs,4*Fs,10*Fs,Fs);
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




