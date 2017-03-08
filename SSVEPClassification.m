 %% SSVEP Classification
% The puspose of this program is to differentiate between different SSVEP
% signals including signal changes. This is accomplished using a variable
% windows. 
    %Clear 
clear;close all;clc;
    %Import Data:
ChannelNames = {['Fp1' 'Fp2' 'Fpz' 'REye']};
% load('mssvep_t2_10_1.mat');
load('mssvep_16.6_3.mat');
remove = 0; % Remove final second of data.
Fs = SamplingRate;
%Import as variables and scale all to one:
removeFromStart = 0;
ch1 = Trial{1}(1+removeFromStart:end-remove,1);
ch2 = Trial{2}(1+removeFromStart:end-remove,1);
ch3 = Trial{3}(1+removeFromStart:end-remove,1);

if size(Trial,2) > 3
    ch4 = Trial{4}(1:end-remove,1);
end
seconds = length(ch1)/Fs
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
% wind = [1024 512 256 128];
[S,wfreqs] = welch_estimator(ch3_f, 250, hann(1024)); 
% S = S(1, :);
figure
plot(wfreqs, S),xlim([1 35]);
showSpect = 1;  
%     showSpect = input('Show Spectrograms?\n');
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
cont = [0];
close all;
minlen = min([ length(ch1) length(ch2) length(ch3) ]);
% winL = [ 250 376 500 626 750 876 1000 ]; %0.5?4s
winL = [250]
newWin = 250;
ii = 1:newWin:(minlen-max(winL));
ops = (length(winL)*length(ii));
xl = [5 25];
wlen = 2^nextpow2(Fs);
h=wlen/4;
win = hamming(wlen, 'periodic'); %CAN BE CONVERTED TO C IF WLEN IS CONSTANT
nfft = 2^nextpow2(wlen+1);
K = sum(hamming(wlen, 'periodic'))/wlen;
if isempty(cont)
    fH = figure(1);
    set(fH, 'Position', [100, 100, 1600, 1000]);
end
loc500 = find(winL==500)-1;
% ----- CLASSIFIER THRESHOLDS ----- %
PeakThresholdPSD = 0.5E-11;
% -- Actual Frequencies: (rounded to 4 decimal places)
f0 = [10.0000 12.5000 15.1515 16.6667];
% Tolerance (f0 +/- f0t)
% Low-res tolerances.
f0t_low = 0.45;
f0r_low = [f0'-f0t_low, f0'+f0t_low];
% High-res tolerances.
f0t = 0.3;
f0r = [f0'-f0t, f0'+f0t];
% Frequency ranges (based on tolerance, for STFT.):
% -- Classifier Possible Outputs:
fp = [10 12 15 16];
tYrejectCount = 1;
for i = 1:length(ii)
    for j = 1:length(winL)
        if isempty(cont)
            fprintf('%d -> %d \n', ii(i),ii(i)+winL(j)-1);
        end
        %TODO: REPLACE CUSTOM FILTER WITH STATIC ONE FOR CONVERSION
        Ch1.Windows{i,j} = customFilt( ch1(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [Ch1.fFFT{j}, Ch1.FFT{i,j}] = get_fft_data(Ch1.Windows{i,j}, Fs);
        [Ch1.MaxFFT{i,j}, Ch1.IndicesMaxFFT{i,j}] = max(Ch1.FFT{i,j});
        
        Ch2.Windows{i,j} = customFilt( ch2(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [~, Ch2.FFT{i,j}] = get_fft_data(Ch2.Windows{i,j}, Fs);
        [Ch2.MaxFFT{i,j}, Ch2.IndicesMaxFFT{i,j}] = max(Ch2.FFT{i,j});
        
        Ch3.Windows{i,j} = customFilt( ch3(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [~, Ch3.FFT{i,j}] = get_fft_data(Ch3.Windows{i,j}, Fs);
        [Ch3.MaxFFT{i,j}, Ch3.IndicesMaxFFT{i,j}] = max(Ch3.FFT{i,j});
            %%% --- POWER SPECTRAL DENSITY EST --- %%%
                % Note: DOES NOT ACCEPT WINDOWS OF ODD LENGTH %
            [Ch1.PSDData{i,j}, Ch1.fPSD{j}] = welch_estimator(Ch1.Windows{i,j}, Fs, hann(winL(j)));
            [Ch2.PSDData{i,j}, Ch2.fPSD{j}] = welch_estimator(Ch2.Windows{i,j}, Fs, hann(winL(j)));
            [Ch3.PSDData{i,j}, Ch3.fPSD{j}] = welch_estimator(Ch3.Windows{i,j}, Fs, hann(winL(j)));
                % FIND MAX VALUE:
            [Ch1.PSDPeak{i,j}, Ch1.PSDi{i,j}] = max(Ch1.PSDData{i,j});
            [Ch2.PSDPeak{i,j}, Ch2.PSDi{i,j}] = max(Ch2.PSDData{i,j});
            [Ch3.PSDPeak{i,j}, Ch3.PSDi{i,j}] = max(Ch3.PSDData{i,j});
                % FIND LOCAL MAX VALUES:
                    % IGNORE EVERYTHING OUTSIDE OF [9 18]Hz
                    % TODO: FIND LOCAL MAX IN FOUR FREQ REGIONS:
            if isempty(cont)
                %%% -------- PLOT PSDs ------------ %%%
                subplot(3,3,[4 6]); hold on;
                plot(Ch1.fPSD{j}, Ch1.PSDData{i,j}),xlim(xl);
                plot(Ch2.fPSD{j}, Ch2.PSDData{i,j}),xlim(xl);
                plot(Ch3.fPSD{j}, Ch3.PSDData{i,j}),xlim(xl);
                % Plot Max Values
                plot(Ch1.fPSD{j}(Ch1.PSDi{i,j}),Ch1.PSDPeak{i,j},'or');
                plot(Ch2.fPSD{j}(Ch2.PSDi{i,j}),Ch2.PSDPeak{i,j},'om');
                plot(Ch3.fPSD{j}(Ch3.PSDi{i,j}),Ch3.PSDPeak{i,j},'og');
                str = [' f = ' num2str(Ch1.fPSD{j}(Ch1.PSDi{i,j}))];
                text(Ch1.fPSD{j}(Ch1.PSDi{i,j}),Ch1.PSDPeak{i,j},str);
                str = [' f = ' num2str(Ch2.fPSD{j}(Ch2.PSDi{i,j}))];
                text(Ch2.fPSD{j}(Ch2.PSDi{i,j}),Ch2.PSDPeak{i,j},str);
                str = [' f = ' num2str(Ch3.fPSD{j}(Ch3.PSDi{i,j}))];
                text(Ch3.fPSD{j}(Ch3.PSDi{i,j}),Ch3.PSDPeak{i,j},str);
                % Plot Local Maxima
                xlabel('Normalized frequency'); %ylabel('PSD [dB]');
                ylabel('Power Spectrum Magnitude');
                title('Power Spectral Test');
                hold off;
                %%% -------- PLOT FFTs ------------ %%%
                subplot(3,3,[1 3]);
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
            end
            %%% --- APPLY STFT --- %%%
            if(winL(j) >= 500)
                %CH1
                [Ch1.sSTFT{i, j-loc500}, Ch1.fSTFT{j-loc500}, Ch1.tSTFT{j-loc500}] = stft( Ch1.Windows{i,j}, h, nfft, Fs );
                Ch1.sSTFT{i, j-loc500} = 20*log10(abs(Ch1.sSTFT{i, j-loc500})/wlen/K + 1e-6); 
                %CH2
                [Ch2.sSTFT{i, j-loc500}, Ch2.fSTFT{j-loc500}, Ch2.tSTFT{j-loc500}] = stft( Ch2.Windows{i,j}, h, nfft, Fs );
                Ch2.sSTFT{i, j-loc500} = 20*log10(abs(Ch2.sSTFT{i, j-loc500})/wlen/K + 1e-6);
                %CH3
                [Ch3.sSTFT{i, j-loc500}, Ch3.fSTFT{j-loc500}, Ch3.tSTFT{j-loc500}] = stft( Ch3.Windows{i,j}, h, nfft, Fs );
                Ch3.sSTFT{i, j-loc500} = 20*log10(abs(Ch3.sSTFT{i, j-loc500})/wlen/K + 1e-6);
                %%%%%%-TODO feature extraction from s f t [5->20Hz]
                if isempty(cont)
                    subplot(3,3,7);
                    imagesc(Ch1.tSTFT{j-loc500},Ch1.fSTFT{j-loc500},Ch1.sSTFT{i, j-loc500}),ylim([7.25 20]);
                    set(gca,'YDir','normal')
                    xlabel('Time, s')
                    ylabel('Frequency, Hz')
                    title('Amplitude spectrogram of the signal')
                    handl = colorbar;
                    ylabel(handl, 'Magnitude, dB')
                    subplot(3,3,8);
                    imagesc(Ch2.tSTFT{j-loc500},Ch2.fSTFT{j-loc500},Ch2.sSTFT{i, j-loc500}),ylim([7.25 20]);
                    set(gca,'YDir','normal')
                    xlabel('Time, s')
                    ylabel('Frequency, Hz')
                    title('Amplitude spectrogram of the signal')
                    handl = colorbar;
                    ylabel(handl, 'Magnitude, dB')
                    subplot(3,3,9);
                    imagesc(Ch3.tSTFT{j-loc500},Ch3.fSTFT{j-loc500},Ch3.sSTFT{i, j-loc500}),ylim([7.25 20]);
                    set(gca,'YDir','normal')
                    xlabel('Time, s')
                    ylabel('Frequency, Hz')
                    title('Amplitude spectrogram of the signal')
                    handl = colorbar;
                    ylabel(handl, 'Magnitude, dB')
                end
            end
        % --- DECISION TREE --- %:
        %1 - Boolean Checks:
        for b = 1:length(f0)
            if winL(j)<=376
                frange = f0r_low(b,:);
                Ch1.B{i,j}.b1(b) = isWithin(Ch1.fPSD{j}(Ch1.PSDi{i,j}),frange);
                Ch2.B{i,j}.b1(b) = isWithin(Ch2.fPSD{j}(Ch1.PSDi{i,j}),frange);
                Ch3.B{i,j}.b1(b) = isWithin(Ch3.fPSD{j}(Ch1.PSDi{i,j}),frange);
                if(Ch1.B{i,j}.b1(b))
                    Ch1.B{i,j}.m1(b) = Ch1.PSDPeak{i,j};
                else
                    Ch1.B{i,j}.m1(b) = 0.0;
                end
                if(Ch2.B{i,j}.b1(b))
                    Ch2.B{i,j}.m1(b) = Ch2.PSDPeak{i,j};
                else
                    Ch2.B{i,j}.m1(b) = 0.0;
                end
                if(Ch3.B{i,j}.b1(b))
                    Ch3.B{i,j}.m1(b) = Ch3.PSDPeak{i,j};
                else
                    Ch3.B{i,j}.m1(b) = 0.0;
                end
            else
                frange = f0r(b,:);
                Ch1.B{i,j}.b1(b) = isWithin(Ch1.fPSD{j}(Ch1.PSDi{i,j}),frange);
                Ch2.B{i,j}.b1(b) = isWithin(Ch2.fPSD{j}(Ch1.PSDi{i,j}),frange);
                Ch3.B{i,j}.b1(b) = isWithin(Ch3.fPSD{j}(Ch1.PSDi{i,j}),frange);
                if(Ch1.B{i,j}.b1(b))
                    Ch1.B{i,j}.m1(b) = Ch1.PSDPeak{i,j};
                else
                    Ch1.B{i,j}.m1(b) = 0.0;
                end
                if(Ch2.B{i,j}.b1(b))
                    Ch2.B{i,j}.m1(b) = Ch2.PSDPeak{i,j};
                else
                    Ch2.B{i,j}.m1(b) = 0.0;
                end
                if(Ch3.B{i,j}.b1(b))
                    Ch3.B{i,j}.m1(b) = Ch3.PSDPeak{i,j};
                else
                    Ch3.B{i,j}.m1(b) = 0.0;
                end
            end
        end
        % --- USER INPUT --- %
        if isempty(cont)
            cont = input('Approve/continue?\n');
            if ~isempty(cont)
                if cont~=0
                    tYreject(tYrejectCount,:) = i;
                    tYrejectCount = tYrejectCount+1;
                    cont = [];
                end
            end
            clf(1);
        end
    end
end

%% Combine Features:
% clearvars -except Ch1 Ch2 Ch3 winL f0
%Preallocate:
tY = zeros(size(Ch1.MaxFFT,1), 1);
for w = 1:length(winL)
    for r = 1:size(Ch1.MaxFFT,1) %Rows
        for c1 = 1
    %         Ch1:
            tXCh1.FFT(r,2*(c1-1)+1) = Ch1.fFFT{c1}(Ch1.IndicesMaxFFT{r,c1});
            tXCh1.FFT(r,2*(c1-1)+2) = Ch1.MaxFFT{r,c1};
            tXCh1.PSD(r,2*(c1-1)+1) = Ch1.fPSD{c1}(Ch1.PSDi{r,c1});
            tXCh1.PSD(r,2*(c1-1)+2) = Ch1.PSDPeak{r,c1};
    %         Ch2:
            tXCh2.FFT(r,2*(c1-1)+1) = Ch1.fFFT{c1}(Ch2.IndicesMaxFFT{r,c1});
            tXCh2.FFT(r,2*(c1-1)+2) = Ch2.MaxFFT{r,c1};
            tXCh2.PSD(r,2*(c1-1)+1) = Ch2.fPSD{c1}(Ch2.PSDi{r,c1});
            tXCh2.PSD(r,2*(c1-1)+2) = Ch2.PSDPeak{r,c1};
    %         Ch3:
            tXCh3.FFT(r,2*(c1-1)+1) = Ch1.fFFT{c1}(Ch3.IndicesMaxFFT{r,c1});
            tXCh3.FFT(r,2*(c1-1)+2) = Ch3.MaxFFT{r,c1};
            tXCh3.PSD(r,2*(c1-1)+1) = Ch3.fPSD{c1}(Ch3.PSDi{r,c1});
            tXCh3.PSD(r,2*(c1-1)+2) = Ch3.PSDPeak{r,c1};
        end
            for c = 1:length(f0)
                tXCh1.M(r, c) = Ch1.B{r,j}.m1(c);
                tXCh2.M(r, c) = Ch2.B{r,j}.m1(c);
                tXCh3.M(r, c) = Ch3.B{r,j}.m1(c);
            end
        tY(r,1) = 16; 
    end
end
tX = [tXCh1.FFT tXCh1.PSD tXCh1.M ...
    tXCh2.FFT tXCh2.PSD tXCh2.M ...
    tXCh3.FFT tXCh3.PSD tXCh3.M ];
tX = [tXCh1 tXCh2 tXCh3];
tXtY = [tX tY];
% clearvars -except tX tY tXtY 

%% Sort out what will be passed to CCA. 
% etx = floor(size(tX,1)/2);
% [Wx, Wy, r] = cca(tX(1:etx,:), tX(etx+1:end,:));



