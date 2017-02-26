%% SSVEP Classification
% The puspose of this program is to differentiate between different SSVEP
% signals including signal changes. This is accomplished using a variable
% windows. 
    %Clear 
clear;close all;clc;
    %Import Data:
ChannelNames = {['Fp1' 'Fp2' 'Fpz' 'REye']};
load('mssvep_12.5_1.mat');
remove = 250; % Remove final second of data.
Fs = SamplingRate;
ch1 = Trial{1}(1:end-remove,1);
ch2 = Trial{2}(1:end-remove,1);
ch3 = Trial{3}(1:end-remove,1);
ch4 = Trial{4}(1:end-remove,1);
flim = [8. 22];
winLim = [6 25];
N = 5;
ch1_f = customFilt(ch1, Fs, flim, N);
ch2_f = customFilt(ch2, Fs, flim, N);
ch3_f = customFilt(ch3, Fs, flim, N);
ch4_f = customFilt(ch4, Fs, flim, N);
[f, P1] = get_fft_data(ch1_f, Fs);
[f2, P2] = get_fft_data(ch2_f, Fs);
[f3, P3] = get_fft_data(ch3_f, Fs);
[f4, P4] = get_fft_data(ch4_f,Fs);
figure(2); hold on;
plot(f,P1,'color','m'),xlim([1 35]);
plot(f2,P2,'color','c'),xlim([1 35]);
plot(f3,P3,'color','r'),xlim([1 35]);
plot(f4,P4,'color','b'),xlim([1 35]);
hold off;
title('FFT(Ch1-4)');
ylabel('|P1(f)|');
xlabel('f (Hz)');
    %Skip Spects. 
   cont = [];
%%-Analysis:
% start with smallest possible window:
minlen = min([length(ch1) length(ch2) length(ch3) length(ch4)]);
winL = [63 125 188 250 313 375 438 500 563 625 688 750 813 875 938 1000 1250 1500 1750 2000]; % Window Lengths (.25s to 4s)
newWin = 250;
ii = 1:newWin:(minlen-max(winL));
ops = (length(winL)*length(ii))
xl = [5 25];

figure(1);
for i = 1:length(ii)
    for j = 1:length(winL)
        Windows.ch1{i,j} = customFilt( ch1(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [FFTData.ch1f{j}, FFTData.ch1fft{i,j}] = get_fft_data(Windows.ch1{i,j}, Fs);
            [FFTData.ch1_M{i,j}, FFTData.ch1_I{i,j}] = max(FFTData.ch1fft{i,j});
        Windows.ch2{i,j} = customFilt( ch2(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [~, FFTData.ch2fft{i,j}] = get_fft_data(Windows.ch2{i,j}, Fs);
            [FFTData.ch2_M{i,j}, FFTData.ch2_I{i,j}] = max(FFTData.ch2fft{i,j});
        Windows.ch3{i,j} = customFilt( ch3(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
        [~, FFTData.ch3fft{i,j}] = get_fft_data(Windows.ch3{i,j}, Fs);
            [FFTData.ch3_M{i,j}, FFTData.ch3_I{i,j}] = max(FFTData.ch3fft{i,j});
%         Windows.ch4{i,j} = customFilt( ch4(ii(i):ii(i)+winL(j)-1), Fs, flim, N);
%         [~, FFTData.ch4fft{i,j}] = get_fft_data(Windows.ch4{i,j}, Fs);
%             [FFTData.ch4_M{i,j}, FFTData.ch4_I{i,j}] = max(FFTData.ch4fft{i,j});
        %%%%% SEE EXAMPLE ON HOW TO APPLY STFT
        % PLOT: 
        if isempty(cont)
            fprintf('%d ? %d \n', ii(i),ii(i)+winL(j)-1);
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
%             plot(FFTData.ch1f{j}, FFTData.ch4fft{i,j}),xlim(xl);
%                 plot(FFTData.ch1f{j}(FFTData.ch4_I{i,j}), FFTData.ch4_M{i,j},'-.b*');
%                 str = [' f = ' num2str(FFTData.ch1f{j}(FFTData.ch4_I{i,j})) '  M = ' num2str( FFTData.ch4_M{i,j} )];
%                 text(FFTData.ch1f{j}(FFTData.ch4_I{i,j}), FFTData.ch4_M{i,j}, str);
            hold off;
            cont = input('continue?\n');
            clf(1);
        end
    end
end