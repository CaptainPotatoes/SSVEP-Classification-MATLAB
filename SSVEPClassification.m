%% SSVEP Classification
% The puspose of this program is to differentiate between different SSVEP
% signals including signal changes. This is accomplished using a variable
% windows. 
    %Clear 
clear;close all;clc;
    %Import Data:
ChannelNames = {['Fp1' 'Fp2' 'Fpz' 'REye']};
load('mssvep_t1_10_1.mat');
remove = 250; % Remove final second of data.
Fs = SamplingRate;
%Import as variables and scale all to one:

ch1 = Trial{1}(remove:end-remove,1);
ch2 = Trial{2}(1:end-remove,1);
ch3 = Trial{3}(1:end-remove,1);
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
 % Skip Spects. 
    cont = [];
 % Spectrogram:
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
   
%% Analysis:
% start with smallest possible window:
minlen = min([ length(ch1) length(ch2) length(ch3) ]);
winL = [63 125 188 250 313 375 438 500 ...
    563 625 688 750 813 875 938 1000 ...
    1250 1500 1750 2000]; % Window Lengths (.25s to 4s)
newWin = 250;
ii = 1:newWin:(minlen-max(winL));
ops = (length(winL)*length(ii));
xl = [5 25];
wlen = 2^nextpow2(Fs);
h=wlen/4;
nfft = 2^nextpow2(wlen+1);
K = sum(hamming(wlen, 'periodic'))/wlen;

figure(1);
for i = 1:length(ii)
    for j = 1:length(winL)
        Windows.ch1{i,j} = customFilt( ch1(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [FFTData.ch1f{j}, FFTData.ch1fft{i,j}] = get_fft_data(Windows.ch1{i,j}, Fs);
            [FFTData.ch1_M{i,j}, FFTData.ch1_I{i,j}] = max(FFTData.ch1fft{i,j});
            if(winL(j) >= 500)
                subplot(2,1,2);
                [STFTData.s{i, j-7}, STFTData.f{i, j-7}, STFTData.t{i, j-7}] = stft( Windows.ch1{i,j}, wlen, h, nfft, Fs );
                STFTData.s{i, j-7} = 20*log10(abs(STFTData.s{i, j-7})/wlen/K + 1e-6); 
                imagesc(STFTData.t{i, j-7}, STFTData.f{i, j-7},STFTData.s{i, j-7}),ylim([7.25 20]);
                %%%%%% TODO window [5 25] 
                set(gca,'YDir','normal')
                set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
                xlabel('Time, s')
                ylabel('Frequency, Hz')
                title('Amplitude spectrogram of the signal')
                handl = colorbar;
                set(handl, 'FontName', 'Times New Roman', 'FontSize', 14)
                ylabel(handl, 'Magnitude, dB')
            end
            
        Windows.ch2{i,j} = customFilt( ch2(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [~, FFTData.ch2fft{i,j}] = get_fft_data(Windows.ch2{i,j}, Fs);
            [FFTData.ch2_M{i,j}, FFTData.ch2_I{i,j}] = max(FFTData.ch2fft{i,j});
        Windows.ch3{i,j} = customFilt( ch3(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [~, FFTData.ch3fft{i,j}] = get_fft_data(Windows.ch3{i,j}, Fs);
            [FFTData.ch3_M{i,j}, FFTData.ch3_I{i,j}] = max(FFTData.ch3fft{i,j});
        %%-%-%% SEE EXAMPLE ON HOW TO APPLY STFT
        % PLOT: 
        if isempty(cont)
            fprintf('%d ? %d \n', ii(i),ii(i)+winL(j)-1);
            subplot(2,1,1);
            hold on;
            plot(FFTData.ch1f{j}, FFTData.ch1fft{i,j}),xlim(xl);
                plot(FFTData.ch1f{j}(FFTData.ch1_I{i,j}), FFTData.ch1_M{i,j},'-.r*');
                str = [' f = ' num2str(FFTData.ch1f{j}(FFTData.ch1_I{i,j})) '  M = ' num2str( FFTData.ch1_M{i,j} )];
                text(FFTData.ch1f{j}(FFTData.ch1_I{i,j}), FFTData.ch1_M{i,j}, str);
            plot(FFTData.ch1f{j}, FFTData.ch2fft{i,j}),xlim(xl);
                 plot(FFTData.ch1f{j}(FFTData.ch2_I{i,j}), FFTData.ch2_M{i,j},'-.m*');
                str = [' f = ' num2str(FFTData.ch1f{j}(FFTData.ch2_I{i,j})) '  M = ' num2str( FFTData.ch2_M{i,j} )];
                text(FFTData.ch1f{j}(FFTData.ch2_I{i,j}), FFTData.ch2_M{i,j}, str);
            plot(FFTData.ch1f{j}, FFTData.ch3fft{i,j}),xlim(xl);
                plot(FFTData.ch1f{j}(FFTData.ch3_I{i,j}), FFTData.ch3_M{i,j},'-.c*');
                str = [' f = ' num2str(FFTData.ch1f{j}(FFTData.ch3_I{i,j})) '  M = ' num2str( FFTData.ch3_M{i,j} )];
                text(FFTData.ch1f{j}(FFTData.ch3_I{i,j}), FFTData.ch3_M{i,j}, str);
            hold off;
            cont = input('continue?\n');
            clf(1);
        end
    end
end