 %% SSVEP Classification
% The puspose of this program is to differentiate between different SSVEP
% signals including signal changes. This is accomplished using a variable
% windows. 
    %Clear 
clear;close all;clc;
    %Import Data:
ChannelNames = {['Fp1' 'Fp2' 'Fpz' 'REye']};
load('mssvep_16.6_3.mat');
remove = 250; % Remove final second of data.
Fs = SamplingRate;
%Import as variables and scale all to one:

ch1 = Trial{1}(1:end-remove,1);
ch2 = Trial{2}(1:end-remove,1);
ch3 = Trial{3}(1:end-remove,1);
% ch3 = ch3(1:2000);
if size(Trial,2) > 3
    ch4 = Trial{4}(1:end-remove,1);
end

flim   = [8.0 18];
winLim = [7.5 20];
N = 5;
    %Filter & Scale everything to '1'
ch1_f = scaleAbs(customFilt(ch1, Fs, flim, N));
ch2_f = scaleAbs(customFilt(ch2, Fs, flim, N));
ch3_f = scaleAbs(customFilt(ch3, Fs, flim, N));
if size(Trial,2) > 3
    ch4_f = scaleAbs(customFilt(ch4, Fs, flim, N));
end
[f,  P1]  = get_fft_data(ch1_f, Fs);
[f2, P2] = get_fft_data(ch2_f, Fs);
[f3, P3] = get_fft_data(ch3_f, Fs);
if size(Trial,2) > 3
    [f4, P4] = get_fft_data(ch4_f, Fs);
end
figure(1); 
hold on;
plot(f,  P1,'color','m'),xlim([1 35]);
plot(f2, P2,'color','c'),xlim([1 35]);
plot(f3, P3,'color','r'),xlim([1 35]);
if size(Trial,2) > 3
    plot(f4, P4,'color','b'),xlim([1 35]);
end
hold off;
title('FFT(Ch1-4)');
ylabel('|P1(f)|');
xlabel('f (Hz)');
wind = [1024 512 256 128];
[S,freqs] = welch_estimator(ch3_f, 250, hann(1024)); 
S = S(1, :);
figure
plot(freqs, S),xlim([1 35]);
cont = [];
 %% Skip Spects.  
%     showSpect = input('Show Spectrograms?\n');
showSpect = 0;
 % Spectrograms:
 if showSpect == 1
    figure;
[~, Fspect, T, P] = spectrogram(ch1_f, 5*Fs,4*Fs,10*Fs,Fs);
imagesc(T, Fspect(Fspect<winLim(2) & Fspect>winLim(1)), 10*log10(P(Fspect<winLim(2) & Fspect>winLim(1),:)));
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
cb = colorbar;
ylabel(cb, 'Power (db)')
colormap(jet)
title('Channel 1', 'FontSize', 14)

    figure;
[~, Fspect, T, P] = spectrogram(ch2_f, 5*Fs,4*Fs,10*Fs,Fs);
imagesc(T, Fspect(Fspect<winLim(2) & Fspect>winLim(1)), 10*log10(P(Fspect<winLim(2) & Fspect>winLim(1),:)));
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
cb = colorbar;
ylabel(cb, 'Power (db)')
colormap(jet)
title('Channel 2', 'FontSize', 14)

    figure;
[~, Fspect, T, P] = spectrogram(ch3_f, 5*Fs,4*Fs,10*Fs,Fs);
imagesc(T, Fspect(Fspect<winLim(2) & Fspect>winLim(1)), 10*log10(P(Fspect<winLim(2) & Fspect>winLim(1),:)));
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel('Time (s)')
cb = colorbar;
ylabel(cb, 'Power (db)')
colormap(jet)
title('Channel 3', 'FontSize', 14)
 end   
%% Analysis (New; Modified):
% start with smallest possible window:
% TODO: EACH winL will result in a feature. USE: FFT, STFT, PSD, and ??
% SEPARATE FEATURES BY 1s, 2s, and 4s elapsed and have separate classifiers
% for each. 
close all;
minlen = min([ length(ch1) length(ch2) length(ch3) ]);
% winL = [63 125 188 250 313 375 438 500 ...
%     563 625 688 750 813 875 938 1000 ...
%     1250 1500 1750 2000]; % Window Lengths (.25s to 4s)
% winL = [64 126 188 250 376 438 500 688 750 938 1000];
winL = [ 126 250 376 500 626 750 876 1000 ];
newWin = 250;
ii = 1:newWin:(minlen-max(winL));
ops = (length(winL)*length(ii));
xl = [5 25];
wlen = 2^nextpow2(Fs);
h=wlen/4;
nfft = 2^nextpow2(wlen+1);
K = sum(hamming(wlen, 'periodic'))/wlen;
if isempty(cont)
    fH = figure(1);
    set(fH, 'Position', [100, 100, 1200, 1000]);
end
loc500 = find(winL==500)-1;
for i = 1:length(ii)
    for j = 1:length(winL)
        Ch1.Windows{i,j} = customFilt( ch1(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [Ch1.fFFT{j}, Ch1.FFT{i,j}] = get_fft_data(Ch1.Windows{i,j}, Fs);
        [Ch1.MaxFFT{i,j}, Ch1.IndicesMaxFFT{i,j}] = max(Ch1.FFT{i,j});
        
        Ch2.Windows{i,j} = customFilt( ch2(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [~, Ch2.FFT{i,j}] = get_fft_data(Ch2.Windows{i,j}, Fs);
        [Ch2.MaxFFT{i,j}, Ch2.IndicesMaxFFT{i,j}] = max(Ch2.FFT{i,j});
        
        Ch3.Windows{i,j} = customFilt( ch3(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [~, Ch3.FFT{i,j}] = get_fft_data(Ch3.Windows{i,j}, Fs);
        [Ch3.MaxFFT{i,j}, Ch3.IndicesMaxFFT{i,j}] = max(Ch3.FFT{i,j});
        
%           --- POWER SPECTRAL DENSITY EST ---
                % Note: DOES NOT ACCEPT WINDOWS OF ODD LENGTH %
            [Ch1.PSDData{i,j}, Ch1.fPSD{j}] = welch_estimator(Ch1.Windows{i,j}, Fs, hann(winL(j)));
            Ch1.PSDData{i,j} = Ch1.PSDData{i,j}(1,:);%3D->2D
            [Ch2.PSDData{i,j}, Ch2.fPSD{j}] = welch_estimator(Ch2.Windows{i,j}, Fs, hann(winL(j)));
            Ch2.PSDData{i,j} = Ch2.PSDData{i,j}(1,:);%3D->2D
            [Ch3.PSDData{i,j}, Ch3.fPSD{j}] = welch_estimator(Ch3.Windows{i,j}, Fs, hann(winL(j)));
            Ch3.PSDData{i,j} = Ch3.PSDData{i,j}(1,:);%3D->2D
            if isempty(cont)
                subplot(3,1,2); hold on;
%                 plot(Ch1.fPSD{i,j}, 10*log10(abs(reshape(PSDData.ch1_PSD{i,j}/2,length(Ch1.fPSD{i,j}),1)))),xlim(xl);
                plot(Ch1.fPSD{j}, Ch1.PSDData{i,j}),xlim(xl);
                plot(Ch2.fPSD{j}, Ch2.PSDData{i,j}),xlim(xl);
                plot(Ch3.fPSD{j}, Ch3.PSDData{i,j}),xlim(xl);
                xlabel('Normalized frequency'); %ylabel('PSD [dB]');
                ylabel('Power Spectrum Magnitude');
                title('Power Spectral Test');
                hold off;
            end
%           --- APPLY STFT ---
            if(winL(j) >= 500)
                subplot(3,1,3);
                %CH1
                [Ch1.sSTFT{i, j-loc500}, Ch1.fSTFT{j-loc500}, Ch1.tSTFT{j-loc500}] = stft( Ch1.Windows{i,j}, wlen, h, nfft, Fs );
                Ch1.sSTFT{i, j-loc500} = 20*log10(abs(Ch1.sSTFT{i, j-loc500})/wlen/K + 1e-6); 
                %CH2
                [Ch2.sSTFT{i, j-loc500}, Ch2.fSTFT{j-loc500}, Ch2.tSTFT{j-loc500}] = stft( Ch2.Windows{i,j}, wlen, h, nfft, Fs );
                Ch2.sSTFT{i, j-loc500} = 20*log10(abs(Ch2.sSTFT{i, j-loc500})/wlen/K + 1e-6);
                %CH3
                [Ch3.sSTFT{i, j-loc500}, Ch3.fSTFT{j-loc500}, Ch3.tSTFT{j-loc500}] = stft( Ch3.Windows{i,j}, wlen, h, nfft, Fs );
                Ch3.sSTFT{i, j-loc500} = 20*log10(abs(Ch3.sSTFT{i, j-loc500})/wlen/K + 1e-6);
                %%%%%%-TODO feature extraction from s f t [5->20Hz]
                if isempty(cont)
                    imagesc(Ch1.tSTFT{j-loc500}, Ch1.fSTFT{j-loc500},Ch1.sSTFT{i, j-loc500}),ylim([7.25 20]);
                    set(gca,'YDir','normal')
                    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
                    xlabel('Time, s')
                    ylabel('Frequency, Hz')
                    title('Amplitude spectrogram of the signal')
                    handl = colorbar;
                    set(handl, 'FontName', 'Times New Roman', 'FontSize', 14)
                    ylabel(handl, 'Magnitude, dB')
%                 imagesc(Ch2.tSTFT{i, j-loc500}, Ch2.fSTFT{i, j-loc500},Ch2.sSTFT{i, j-loc500}),ylim([7.25 20]);
%                 imagesc(Ch3.tSTFT{i, j-loc500}, Ch3.fSTFT{i, j-loc500},Ch3.sSTFT{i, j-loc500}),ylim([7.25 20]);
                end
            end
        % PLOT: 
        if isempty(cont)
            fprintf('%d -> %d \n', ii(i),ii(i)+winL(j)-1);
            subplot(3,1,1);
            hold on;
            plot(Ch1.fFFT{j}, Ch1.FFT{i,j}),xlim(xl);
                plot(Ch1.fFFT{j}(Ch1.IndicesMaxFFT{i,j}), Ch1.MaxFFT{i,j},'-.r*');
                str = [' f = ' num2str(Ch1.fFFT{j}(Ch1.IndicesMaxFFT{i,j})) '  M = ' num2str( Ch1.MaxFFT{i,j} )];
                text(Ch1.fFFT{j}(Ch1.IndicesMaxFFT{i,j}), Ch1.MaxFFT{i,j}, str);
            plot(Ch1.fFFT{j}, Ch2.FFT{i,j}),xlim(xl);
                plot(Ch1.fFFT{j}(Ch2.IndicesMaxFFT{i,j}), Ch2.MaxFFT{i,j},'-.m*');
                str = [' f = ' num2str(Ch1.fFFT{j}(Ch2.IndicesMaxFFT{i,j})) '  M = ' num2str( Ch2.MaxFFT{i,j} )];
                text(Ch1.fFFT{j}(Ch2.IndicesMaxFFT{i,j}), Ch2.MaxFFT{i,j}, str);
            plot(Ch1.fFFT{j}, Ch3.FFT{i,j}),xlim(xl);
                plot(Ch1.fFFT{j}(Ch3.IndicesMaxFFT{i,j}), Ch3.MaxFFT{i,j},'-.c*');
                str = [' f = ' num2str(Ch1.fFFT{j}(Ch3.IndicesMaxFFT{i,j})) '  M = ' num2str( Ch3.MaxFFT{i,j} )];
                text(Ch1.fFFT{j}(Ch3.IndicesMaxFFT{i,j}), Ch3.MaxFFT{i,j}, str);
            title('FFT (Ch 1-3): With Peaks');
            ylabel('|P1(f)|');
            xlabel('f (Hz)');
            hold off;
            cont = input('continue?\n');
            clf(1);
        end
    end
end
clearvars -except Ch1 Ch2 Ch3

