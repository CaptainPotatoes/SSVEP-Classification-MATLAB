function varargout = SSVEPGUI(varargin)
% SSVEPGUI MATLAB code for SSVEPGUI.fig
%      SSVEPGUI, by itself, creates a new SSVEPGUI or raises the existing
%      singleton*.
%
%      H = SSVEPGUI returns the handle to a new SSVEPGUI or the handle to
%      the existing singleton*.
%
%      SSVEPGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SSVEPGUI.M with the given input arguments.
%
%      SSVEPGUI('Property','Value',...) creates a new SSVEPGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the SSVEPGUI before SSVEPGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SSVEPGUI_OpeningFcn via varargin.
%
%      *See SSVEPGUI Options on GUIDE's Tools menu.  Choose "SSVEPGUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help SSVEPGUI
% Last Modified by GUIDE v2.5 27-Feb-2017 13:16:51
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SSVEPGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SSVEPGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before SSVEPGUI is made visible.
function SSVEPGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SSVEPGUI (see VARARGIN)

% global countFar countMiddle countClose
% Choose default command line output for SSVEPGUI
handles.output = hObject;
global trainingData
trainingData = cell(1,2);
trainingData{1} = [0,0];
trainingData{2} = [0,0];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SSVEPGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = SSVEPGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global myDevice;
global deviceName;
global which_pc;
count = 0;

if get(hObject,'Value') && count == 0
which_pc = input('WHICH PC? : 1=home, 0=work \n');
% load the BioRadio API using a MATLAB's .NET interface
if which_pc ==0
    [ deviceManager , flag ] = load_API('C:\Users\mahmoodms\Dropbox\Public\_VCU\Yeo Lab\_SSVEP\_MATLAB-SSVEP-Classification\BioRadioSDK.dll');
elseif which_pc ==1
    [ deviceManager , flag ] = load_API('C:\Users\Musa Mahmood\Dropbox\Public\_VCU\Yeo Lab\_SSVEP\_MATLAB-SSVEP-Classification\BioRadioSDK.dll');
end
% input = full path to api dll file
% outputs = deviceManager object, success flag

if ~flag % if API not successfully loaded, do not continue
    return
end

% search for available sensors and select one
%
[ deviceName , macID , ok ] = BioRadio_Find( deviceManager );
% input = deviceManager object
% outputs = device name, macid, and flag if selection was canceled out
%

if ~ok %if no sensors selected, do not continue
    errordlg('Please select a BioRadio.')
    return    
% initialize BioRadio object
end


[ myDevice, flag ] = BioRadio_Connect ( deviceManager , macID , deviceName );
% input = deviceManager object, 64-bit mac address of BioRadio, and name of
% BioRadio
% outputs = BioRadio object, success flag for connection
%global myDevice;
if ~flag %if connection failed, do not continue
    return
end
count = 1;

else 
    BioRadio_Disconnect( myDevice )
    count = 0;
    
end
% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global myDevice Idx 
global trainingData totalCount which_pc

totalCount = cell(2);
totalCount{1} = 0;
totalCount{2} = 0;
BioRadio_Name = 'EEG-SSVEP';
numEnabledBPChannels = double(myDevice.BioPotentialSignals.Count);

if numEnabledBPChannels == 0
    myDevice.Disconnect;
    RawBioRadioData = [];
    errordlg('No BioPotential Channels Programmed. Return to BioCapture to Configure.')
    return
end

sampleRate_BP = double(myDevice.BioPotentialSignals.SamplesPerSecond);
%Preallocating and setting up which area in the SSVEPGUI the plot will go into
    % First two are raw data
    % Next two are data analysis features. 
numAxes = 14; 
axis_handles = zeros(1,numAxes);
for ch = 1:numAxes
    axis_handles(ch) = handles.(['axes',num2str(ch)]);
end

%Preallocating BPSignals
BioPotentialSignals = cell(1,numEnabledBPChannels);
filteredEEG{ch} = cell(1,numEnabledBPChannels);
Idx = cell(1,numEnabledBPChannels);
if get(hObject,'Value') == 1
myDevice.StartAcquisition;
end
plotWindow = 5;
plotGain_BP = 1;
fft_len = plotWindow*sampleRate_BP;
% fft_len = 2^(nextpow2(plotWindow*sampleRate_BP)); 
% USE WITH fft(X,fft_len_pow2)
f0 = [8.75 18.5];
xl = [5 20];
N = 5;
spect_1 = handles.axes5;
spect_2 = handles.axes6;
t = cell(numEnabledBPChannels, 1);
fp1_data_unfilt = zeros(plotWindow*sampleRate_BP,1);
fp2_data_unfilt = zeros(plotWindow*sampleRate_BP,1);
eog3_data_unfilt = zeros(plotWindow*sampleRate_BP,1); 
eog4_data_unfilt = zeros(plotWindow*sampleRate_BP,1);
%LOAD tXtY:: (NOW EMBEDDED IN CODE AS STATIC VAR)
% load('allEOGtD.mat');
% SSVEP Training Data:
% load('tXtY_SSVEP.mat');
ln = 250-1;
for i = 1:numEnabledBPChannels
    cIdx{i} = 1;
end
Y = cell(5,1);
OUTPUT = cell(1,1);
op=1;
while get(hObject,'Value') == 1
    pause(0.05)
    for ch = 1:numEnabledBPChannels
            BioPotentialSignals{ch} = [BioPotentialSignals{ch}; myDevice.BioPotentialSignals.Item(ch-1).GetScaledValueArray.double'];
            Idx{ch} = 1:length(BioPotentialSignals{ch});            
            %Plot the Axes in the SSVEPGUI
            if length(BioPotentialSignals{ch}) <= plotWindow*sampleRate_BP
                t{ch} = (0:(length(BioPotentialSignals{ch})-1))*(1/sampleRate_BP);
                if ch==1 || ch==2
                    if length(BioPotentialSignals{ch}) >= 31
                        filteredEEG{ch} = eeg_h_custom(plotGain_BP*BioPotentialSignals{ch}, sampleRate_BP, f0, N);
                        plot(axis_handles(ch),t{ch},filteredEEG{ch});
                    end
                    set(handles.(['axes',num2str(ch)]),'XLim',[0 plotWindow]);
                    set(get(handles.(['axes',num2str(ch)]), 'XLabel'), 'String', 'Time(s)')
                    set(get(handles.(['axes',num2str(ch)]), 'YLabel'), 'String',  'mV')
                end
                if ch==1
                    set(get(handles.(['axes',num2str(ch)]), 'Title'), 'String', 'Channel 1 Data')
                    if length(BioPotentialSignals{ch})>sampleRate_BP
                        fp1_data_filtered = eeg_h_custom(BioPotentialSignals{ch}, sampleRate_BP, f0, N);
                        [f, P1] = get_fft_data(fp1_data_filtered, sampleRate_BP);
                        plot(axis_handles(3),f,P1);
                        set(handles.(['axes',num2str(3)]),'XLim',xl);
                        set(get(handles.(['axes',num2str(3)]), 'XLabel'), 'String', 'f (Hz)')
                        set(get(handles.(['axes',num2str(3)]), 'YLabel'), 'String', '|P1(f)|')
                        set(get(handles.(['axes',num2str(3)]), 'Title'), 'String', 'FFT(Fp1)')
                    end
                elseif ch==2
                    set(get(handles.(['axes',num2str(ch)]),'Title'),'String','Channel 2 Data')
                    if length(BioPotentialSignals{ch})>sampleRate_BP
                        fp2_data_filtered = eeg_h_custom(BioPotentialSignals{ch}, sampleRate_BP, f0, N);
                        [f, P1] = get_fft_data(fp2_data_filtered, sampleRate_BP);
                        plot(axis_handles(4),f,P1);
                        set(handles.(['axes',num2str(4)]),'XLim',xl);
                        set(get(handles.(['axes',num2str(4)]), 'XLabel'), 'String', 'f (Hz)')
                        set(get(handles.(['axes',num2str(4)]), 'YLabel'), 'String', '|P1(f)|')
                        set(get(handles.(['axes',num2str(4)]), 'Title'), 'String', 'FFT(Fp2)')
                    end
                end
            else %once plot window is exceeded:
                if ch==1
                     t{1} = ((length(BioPotentialSignals{ch})-(plotWindow*sampleRate_BP-1)):length(BioPotentialSignals{ch}))*(1/sampleRate_BP);
                     filteredEEG{ch} = eeg_h_custom(plotGain_BP*BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end), sampleRate_BP, f0, N);
                     plot(axis_handles(ch),t{1},filteredEEG{ch})
                     set(get(handles.(['axes',num2str(ch)]), 'Title'), 'String', 'Channel 1 Filtered')
                elseif ch==2
                     t{2} = ((length(BioPotentialSignals{ch})-(plotWindow*sampleRate_BP-1)):length(BioPotentialSignals{ch}))*(1/sampleRate_BP);
                     filteredEEG{ch} = eeg_h_custom(plotGain_BP*BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end), sampleRate_BP, f0, N);
                     plot(axis_handles(ch),t{2},filteredEEG{ch})
                     set(get(handles.(['axes',num2str(ch)]), 'Title'), 'String', 'Channel 2 Filtered')
                elseif ch==3
                     t{3} = ((length(BioPotentialSignals{ch})-(plotWindow*sampleRate_BP-1)):length(BioPotentialSignals{ch}))*(1/sampleRate_BP);
                elseif ch==4
                     t{4} = ((length(BioPotentialSignals{ch})-(plotWindow*sampleRate_BP-1)):length(BioPotentialSignals{ch}))*(1/sampleRate_BP);
                end
                if ch == 1 || ch == 2
                    set(handles.(['axes',num2str(ch)]),'XLim',[t{1}(end)-plotWindow t{1}(end)]);
                    set(handles.(['axes',num2str(ch)]),'YLim',[-2E-4 2E-4]);
                    set(get(handles.(['axes',num2str(ch)]), 'XLabel'), 'String', 'Time(s)')
                    set(get(handles.(['axes',num2str(ch)]), 'YLabel'), 'String',  'mV')
                end
                if ch==1
                    % Window and Filter
                    fp1_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end);
                    fp1_data_filtered = eeg_h_custom(fp1_data_unfilt, sampleRate_BP, f0, N);
                    [f, P1] = get_fft_data(fp1_data_filtered, sampleRate_BP);
                    plot(axis_handles(3),f,P1);
                    set(handles.(['axes',num2str(3)]),'XLim',xl);
                    set(get(handles.(['axes',num2str(3)]), 'XLabel'), 'String', 'f (Hz)')
                    set(get(handles.(['axes',num2str(3)]), 'YLabel'), 'String', '|P1(f)|')
                    set(get(handles.(['axes',num2str(3)]), 'Title'), 'String', 'FFT(Fp1)')
                    % Spectrogram (STFT):
                    if which_pc == 0
%                         [S, Fspect, T, P] = spectrogram(fp1_data_filtered, 5*sampleRate_BP,4*sampleRate_BP,10*sampleRate_BP,sampleRate_BP);
%                         imagesc(spect_1, T, Fspect(Fspect<30 & Fspect>1), 10*log10(P(Fspect<30 & Fspect>1,:)));
%                         set(spect_1,'YDir','normal')
%                         cb = colorbar(spect_1);
%                         ylabel(cb, 'Power (db)')
%                         colormap(spect_1,jet)
%                         set(get(handles.(['axes',num2str(5)]), 'XLabel'), 'String', 'Time (s)')
%                         set(get(handles.(['axes',num2str(5)]), 'YLabel'), 'String', 'Frequency (Hz)')
%                         set(get(handles.(['axes',num2str(5)]), 'Title'), 'String', 'Spectrogram (Fp1)')
                    end
                    % Pwelch
                    [Pxx, F] = pwelch(fp1_data_filtered,[],[],250);
                    plot(axis_handles(5), Pxx); 
                    set(handles.(['axes',num2str(5)]),'XLim',xl);
                    set(get(handles.(['axes',num2str(5)]), 'XLabel'), 'String', 'Frequency (Hz)')
                    set(get(handles.(['axes',num2str(5)]), 'YLabel'), 'String', 'Power (dB)')
                    set(get(handles.(['axes',num2str(5)]), 'Title'), 'String', 'Pwelch (Fp1)')
                    %EOG filt:
%                     eog1_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end);
                    eog_data_1 = eog_h_fcn(fp1_data_unfilt, sampleRate_BP);
                    plot(axis_handles(9),t{ch},eog_data_1);
                    set(handles.(['axes',num2str(9)]),'XLim',[t{1}(end)-plotWindow t{1}(end)]);
                    set(handles.(['axes',num2str(9)]),'YLim',[-4E-4 4E-4]);
                    set(get(handles.(['axes',num2str(9)]), 'XLabel'), 'String', 'Time (s)')
                    set(get(handles.(['axes',num2str(9)]), 'YLabel'), 'String',  'mV')
                    set(get(handles.(['axes',num2str(9)]), 'Title'), 'String', ['EOG Filter Ch' num2str(ch)])
                elseif ch==2
                    fp2_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end);
                    fp2_data_filtered = eeg_h_custom(fp2_data_unfilt, sampleRate_BP, f0, N);
                    [f, P1] = get_fft_data(fp2_data_filtered, sampleRate_BP);
                    plot(axis_handles(4),f,P1);
                    set(handles.(['axes',num2str(4)]),'XLim',xl);
                    set(get(handles.(['axes',num2str(4)]), 'XLabel'), 'String', 'f (Hz)')
                    set(get(handles.(['axes',num2str(4)]), 'YLabel'), 'String', '|P2(f)|')
                    set(get(handles.(['axes',num2str(4)]), 'Title'), 'String', 'FFT(Fp2)')
                    % Spect:
                    if which_pc == 0
%                         [S, Fspect, T, P] = spectrogram(fp2_data_filtered, 5*sampleRate_BP,4*sampleRate_BP,10*sampleRate_BP,sampleRate_BP);
%                         imagesc(spect_2, T, Fspect(Fspect<30 & Fspect>1), 10*log10(P(Fspect<30 & Fspect>1,:)));
%                         set(spect_2,'YDir','normal') %gca
%                         cb2 = colorbar(spect_2);
%                         ylabel(cb2, 'Power (db)')
%                         colormap(spect_2,jet)
%                         set(get(handles.(['axes',num2str(6)]), 'XLabel'), 'String', 'Time (s)')
%                         set(get(handles.(['axes',num2str(6)]), 'YLabel'), 'String', 'Frequency (Hz)')
%                         set(get(handles.(['axes',num2str(6)]), 'Title'), 'String', 'Spectrogram (Fp2)')
                    end
                    % P welch
                    [Pxx, F] = pwelch(fp2_data_filtered,[],[],250);
                    plot(axis_handles(6), Pxx); 
                    set(handles.(['axes',num2str(6)]),'XLim',xl);
                    set(get(handles.(['axes',num2str(6)]), 'XLabel'), 'String', 'Frequency (Hz)')
                    set(get(handles.(['axes',num2str(6)]), 'YLabel'), 'String', 'Power (dB)')
                    set(get(handles.(['axes',num2str(6)]), 'Title'), 'String', 'Pwelch (Fp2)')
                    %EOG Filt
%                     eog2_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end);
                    eog_data_1 = eog_h_fcn(fp2_data_unfilt, sampleRate_BP);
                    plot(axis_handles(10),t{ch},eog_data_1);
                    set(handles.(['axes',num2str(10)]),'XLim',[t{1}(end)-plotWindow t{1}(end)]);
                    set(handles.(['axes',num2str(10)]),'YLim',[-4E-4 4E-4]);
                    set(get(handles.(['axes',num2str(10)]), 'XLabel'), 'String', 'Time (s)')
                    set(get(handles.(['axes',num2str(10)]), 'YLabel'), 'String',  'mV')
                    set(get(handles.(['axes',num2str(10)]), 'Title'), 'String', ['EOG Filter Ch' num2str(ch)])
                end
                if ch==3
                    eog3_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end);
                    eog_data_1 = eog_h_fcn(eog3_data_unfilt, sampleRate_BP);
                    plot(axis_handles(7),t{ch},eog_data_1);
                    set(handles.(['axes',num2str(7)]),'XLim',[t{1}(end)-plotWindow t{1}(end)]);
                    set(handles.(['axes',num2str(7)]),'YLim',[-4E-4 4E-4]);
                    set(get(handles.(['axes',num2str(7)]), 'XLabel'), 'String', 'Time (s)')
                    set(get(handles.(['axes',num2str(7)]), 'YLabel'), 'String',  'mV')
                    set(get(handles.(['axes',num2str(7)]), 'Title'), 'String', ['EOG Filter Ch' num2str(ch)])
                    %%FFT CH3:
                    CH3_data_filtered = eeg_h_custom(eog3_data_unfilt, sampleRate_BP, f0, N);
                    [f, P1] = get_fft_data(CH3_data_filtered, sampleRate_BP);
                    plot(axis_handles(11),f,P1);
                    set(handles.(['axes',num2str(11)]),'XLim',xl);
                    set(get(handles.(['axes',num2str(11)]), 'XLabel'), 'String', 'f (Hz)')
                    set(get(handles.(['axes',num2str(11)]), 'YLabel'), 'String', '|P1(f)|')
                    set(get(handles.(['axes',num2str(11)]), 'Title'), 'String', 'FFT(Ch3, Fpz)')
                    %PSD CH3:
                    [Pxx, F] = pwelch(CH3_data_filtered,[],[],250);
                    plot(axis_handles(12), Pxx); 
                    set(handles.(['axes',num2str(12)]),'XLim',xl);
                    set(get(handles.(['axes',num2str(12)]), 'XLabel'), 'String', 'Frequency (Hz)')
                    set(get(handles.(['axes',num2str(12)]), 'YLabel'), 'String', 'Power (dB)')
                    set(get(handles.(['axes',num2str(12)]), 'Title'), 'String', 'Pwelch (Ch3)')
                elseif ch==4
                    eog4_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end);
                    eog_data_2 = eog_h_fcn(eog4_data_unfilt, sampleRate_BP);
                    plot(axis_handles(8),t{ch},eog_data_2);
                    set(handles.(['axes',num2str(8)]),'XLim',[t{2}(end)-plotWindow t{2}(end)]);
                    set(handles.(['axes',num2str(8)]),'YLim',[-4E-4 4E-4]);
                    set(get(handles.(['axes',num2str(8)]), 'XLabel'), 'String', 'Time (s)')
                    set(get(handles.(['axes',num2str(8)]), 'YLabel'), 'String',  'mV')
                    set(get(handles.(['axes',num2str(8)]), 'Title'), 'String', ['EOG Filter Ch' num2str(ch)])
                    %ffT
                    CH4_data_filtered = eeg_h_custom(eog4_data_unfilt, sampleRate_BP, f0, N);
                    [f, P1] = get_fft_data(CH4_data_filtered, sampleRate_BP);
                    plot(axis_handles(13),f,P1);
                    set(handles.(['axes',num2str(13)]),'XLim',xl);
                    set(get(handles.(['axes',num2str(13)]), 'XLabel'), 'String', 'f (Hz)')
                    set(get(handles.(['axes',num2str(13)]), 'YLabel'), 'String', '|P1(f)|')
                    set(get(handles.(['axes',num2str(13)]), 'Title'), 'String', 'FFT(Ch4, Fpz)')
                    %PSD CH4:
                    [Pxx, F] = pwelch(CH4_data_filtered,[],[],250); 
                    plot(axis_handles(14), Pxx); 
                    set(handles.(['axes',num2str(14)]),'XLim',xl);
                    set(get(handles.(['axes',num2str(14)]), 'XLabel'), 'String', 'Frequency (Hz)')
                    set(get(handles.(['axes',num2str(14)]), 'YLabel'), 'String', 'Power (dB)')
                    set(get(handles.(['axes',num2str(14)]), 'Title'), 'String', 'Pwelch (Ch4)')
                end
%                 CLASSIFY EOG DATA:
                   %take each ch unfilt (5s) & cut into 5 pieces, & pass
                   %thru fullHybridClassifier: Put in columns:
                for i = 1:numEnabledBPChannels
                    b1(i) = length(BioPotentialSignals{i}) > cIdx{i}+125;
                end
                
                if sum(b1) == numEnabledBPChannels
                    for i = 1:numEnabledBPChannels
                        cIdx{i} = length(BioPotentialSignals{i});       %Assign new current indx
                        W{i} = BioPotentialSignals{i}(end-ln:end);
                    end

                    if ln<499%ln==249
                        Y{1} = fullHybridClassifier2(W{1}, W{2}, W{3}, W{4}, sampleRate_BP)';
                        disp(Y{1})
                        ln = ln+250;
                    elseif ln >= 499 && ln <749 %ln==499
                        Y{2} = fullHybridClassifier2(W{1}, W{2}, W{3}, W{4}, sampleRate_BP)';
                        disp(Y{2})
                        ln = ln+250;
                    elseif ln>=740 && ln<999%ln == 749
                        Y{3} = fullHybridClassifier2(W{1}, W{2}, W{3}, W{4}, sampleRate_BP)';
                        disp(Y{3})
                        ln = ln+250;
                    elseif ln >= 999 && ln <1249
                        Y{4} = fullHybridClassifier2(W{1}, W{2}, W{3}, W{4}, sampleRate_BP)';
                        disp(Y{4})
                        ln = ln+250;
                    elseif ln >= 1249
                        Y{5} = fullHybridClassifier2(W{1}, W{2}, W{3}, W{4}, sampleRate_BP)';
                        disp(Y{5})
                        
                        if (length(Y{5})==5)
                            if(Y{5}(2) == Y{5}(3)) && (Y{5}(5) == 1)
                                OUTPUT{1}(op) = Y{5}(2);
                                ln=249;
                            else
                                OUTPUT{1}(op) = 0;
                                ln=ln+250;
                            end
                        else
                            OUTPUT{1}(op) = Y{5}(1);
                            ln = ln+250;
                        end
                        fprintf('\n >>>>OUTPUT = %d\n',OUTPUT{1}(op));
                        op=op+1; 
                    end
                    for y = 1:5
                        if ~isempty(Y{y})
                            if (Y{y}(1) == 1)
                                %Double Blink Detected: Reset Command (Stop and
                                %reset):
                                Y = cell(5,1); %Reset cell
                                ln = 249;
                                fprintf('\n >>>>OUTPUT = DOUBLE BLINK :: RESET COMMAND \n \n');
                            end
                        end
                    end
                end
            end

            %% Todo: FFT and plot on axes #3
%             if length(BioPotentialSignals{ch}) > 500
%             plot(axis_handles(ch+1),lags{ch}(size(BufferFilt{1},2)-1,:)/sampleRate_BP,BPAutocorrelation{ch}(size(BufferFilt{1},2)-1,:))
%             end
            
    end     %/for ch = 1:numEnabledBPChannels
end     %/while connected==1

if get(hObject,'Value') == 0
    myDevice.StopAcquisition;
    
    Trial = cell(1,numEnabledBPChannels);
    for i=1:numEnabledBPChannels
        Trial{1,i} = BioPotentialSignals{i};
    end
    assignin('base','NumberOfChannels',numEnabledBPChannels);
    assignin('base','Trial',Trial)
    SamplingRate = sampleRate_BP;
    assignin('base','SamplingRate',SamplingRate);
    %% Change into button function.
    OP_1 = OUTPUT{1};
    H_Notes = [handles.edit2 handles.edit_ChannelLocations, handles.editImpedanceValues, ...
        handles.editElectrodeType, handles.editStimulusSource, handles.editMiscInfo];
    RecordingNotes = cell(length(H_Notes)+1, 2);
    RecordingNotes{1,1} = 'Filename (.mat)'; RecordingNotes{2,1} = 'Channel Locations';
    RecordingNotes{3,1} = 'Electrode Impedance'; RecordingNotes{4,1} = 'Electrode Type';
    RecordingNotes{5,1} = 'Source of Stimulus'; RecordingNotes{6,1} = 'Misc Notes';
    for i=1:length(H_Notes)
        RecordingNotes{i,2} = get(H_Notes(i),'String');
        if isempty(RecordingNotes{i,2})
            RecordingNotes{i,2} = 'No Notes Recorded for This Session'; 
        end
    end
    RecordingNotes{length(H_Notes)+1,1} = 'Sampling Rate';
    RecordingNotes{length(H_Notes)+1,2} = num2str(SamplingRate);
    assignin('base','RecordingNotes',RecordingNotes);
    assignin('base','OUTPUT',OP_1);
    %% Auto-save variables
    filename = get(handles.edit2,'String');
    if isempty(filename)
        c=clock;
        filename = ['RawBioRadioData_',num2str(c(2)),'-',num2str(c(3)),'-',num2str(c(1)),'_',num2str(c(4)),'.',num2str(c(5)),',',num2str(c(6)),'.mat'];
    end
    save(filename,'Trial','SamplingRate','RecordingNotes','trainingData','OUTPUT');
end

% --- Executes on button press in pushbuttonDoubleBlink.
function pushbuttonDoubleBlink_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDoubleBlink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx doubleBlinkIdx doubleBlinkCount trainingData totalCount
if doubleBlinkCount==1
    doubleBlinkIdx = zeros(1,1);
    
    doubleBlinkIdx(1,1) = size(Idx{1,1},2)
    
    doubleBlinkCount = doubleBlinkCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = doubleBlinkIdx; %index
    trainingData{1}(totalCount{1},2) = 1; %class
    
    assignin('base','tD',trainingData);
else
    doubleBlinkIdx(1,1) = size(Idx{1,1},2)
    
    doubleBlinkCount = doubleBlinkCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = doubleBlinkIdx; %index
    trainingData{1}(totalCount{1},2) = 1; %class
    
    assignin('base','tD',trainingData);
end

% --- Executes on button press in pushbuttonSingleBlink.
function pushbuttonSingleBlink_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSingleBlink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx singleBlinkIdx singleBlinkCount trainingData totalCount
if singleBlinkCount==1
    singleBlinkIdx = zeros(1,1);
    
    singleBlinkIdx(1,1) = size(Idx{1,1},2)
    
    singleBlinkCount = singleBlinkCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = singleBlinkIdx; %index
    trainingData{1}(totalCount{1},2) = -1; %class
    
    assignin('base','tD',trainingData);
else
    singleBlinkIdx(1,1) = size(Idx{1,1},2)
    
    singleBlinkCount = singleBlinkCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = singleBlinkIdx; %index
    trainingData{1}(totalCount{1},2) = -1; %class
   
    assignin('base','tD',trainingData);
end

% --- Executes on button press in pushbuttonUp.
function pushbuttonUp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx eyeMoveUpIdx eyeMoveCount trainingData totalCount
if eyeMoveCount==1
    eyeMoveUpIdx = zeros(1,1);
    
    eyeMoveUpIdx(1,1) = size(Idx{1,1},2)
    
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = eyeMoveUpIdx; %index
    trainingData{1}(totalCount{1},2) = 2; %class
    assignin('base','tD',trainingData);
else
    eyeMoveUpIdx(1,1) = size(Idx{1,1},2)
    
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = eyeMoveUpIdx; %index
    trainingData{1}(totalCount{1},2) = 2; %class
    
    assignin('base','tD',trainingData);
end

% --- Executes on button press in pushbuttonDown.
function pushbuttonDown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx eyeMoveDownIdx eyeMoveCount trainingData totalCount
if eyeMoveCount==1
    eyeMoveDownIdx = zeros(1,1);
    
    eyeMoveDownIdx(1,1) = size(Idx{1,1},2)
    
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = eyeMoveDownIdx; %index
    trainingData{1}(totalCount{1},2) = 3; %class
    assignin('base','tD',trainingData);
else
    eyeMoveDownIdx(1,1) = size(Idx{1,1},2)
    
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = eyeMoveDownIdx; %index
    trainingData{1}(totalCount{1},2) = 3; %class
    
    assignin('base','tD',trainingData);
end

% --- Executes on button press in pushbuttonLeft.
function pushbuttonLeft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx eyeMoveLeftIdx eyeMoveCount trainingData totalCount
if eyeMoveCount==1
    eyeMoveLeftIdx = zeros(1,1);
    
    eyeMoveLeftIdx(1,1) = size(Idx{1,1},2)
    
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = eyeMoveLeftIdx; %index
    trainingData{1}(totalCount{1},2) = 4; %class
    assignin('base','tD',trainingData);
else
    eyeMoveLeftIdx(1,1) = size(Idx{1,1},2)
    
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = eyeMoveLeftIdx; %index
    trainingData{1}(totalCount{1},2) = 4; %class
    
    assignin('base','tD',trainingData);
end

% --- Executes on button press in pushbuttonRight.
function pushbuttonRight_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx eyeMoveRightIdx eyeMoveCount trainingData totalCount
if eyeMoveCount==1
    eyeMoveRightIdx = zeros(1,1);
    
    eyeMoveRightIdx(1,1) = size(Idx{1,1},2)
    
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = eyeMoveRightIdx; %index
    trainingData{1}(totalCount{1},2) = 5; %class
    assignin('base','tD',trainingData);
else
    eyeMoveRightIdx(1,1) = size(Idx{1,1},2)
    
    eyeMoveCount = eyeMoveCount+1;
    totalCount{1} = totalCount{1}+1;
    
    trainingData{1}(totalCount{1},1) = eyeMoveRightIdx; %index
    trainingData{1}(totalCount{1},2) = 5; %class
    
    assignin('base','tD',trainingData);
end


% --- Executes on button press in pushbutton11. (10Hz)
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx freq10Idx freqCount trainingData totalCount
if freqCount==1
    freq10Idx = zeros(1,1);
    freq10Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
    
    trainingData{2}(totalCount{2},1) = freq10Idx;
    trainingData{2}(totalCount{2},2) = 1; %class 
    assignin('base','tD',trainingData);
else
    freq10Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
    
    trainingData{2}(totalCount{2},1) = freq10Idx;
    trainingData{2}(totalCount{2},2) = 1; %class 
    assignin('base','tD',trainingData);
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx freq12Idx freqCount trainingData totalCount
if freqCount==1
    freq12Idx = zeros(1,1);
    freq12Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
    
    trainingData{2}(totalCount{2},1) = freq12Idx;
    trainingData{2}(totalCount{2},2) = 2; %class 
    assignin('base','tD',trainingData);
else
    freq12Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
    
    trainingData{2}(totalCount{2},1) = freq12Idx;
    trainingData{2}(totalCount{2},2) = 2; %class 
    assignin('base','tD',trainingData);
end

% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx freq15Idx freqCount trainingData totalCount
if freqCount==1
    freq15Idx = zeros(1,1);
    freq15Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
    
    trainingData{2}(totalCount{2},1) = freq15Idx;
    trainingData{2}(totalCount{2},2) = 3; %class 
    assignin('base','tD',trainingData);
else
    freq15Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
    
    trainingData{2}(totalCount{2},1) = freq15Idx;
    trainingData{2}(totalCount{2},2) = 3; %class 
    assignin('base','tD',trainingData);
end
% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Idx freq16Idx freqCount trainingData totalCount

if freqCount==1
    freq16Idx = zeros(1,1);
    freq16Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
    
    trainingData{2}(totalCount{2},1) = freq16Idx;
    trainingData{2}(totalCount{2},2) = 3; %class 
    assignin('base','tD',trainingData);
else
    freq16Idx(1,1) = size(Idx{1,1},2)
    freqCount = freqCount+1;
    totalCount{2} = totalCount{2}+1;
    
    trainingData{2}(totalCount{2},1) = freq16Idx;
    trainingData{2}(totalCount{2},2) = 3; %class 
    assignin('base','tD',trainingData);
end


%------------------------ %%%%%%%%%%%% ?IGNORE? %%%%%%%%%%% ---------------------------%

function edit_ChannelLocations_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ChannelLocations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ChannelLocations as text
%        str2double(get(hObject,'String')) returns contents of edit_ChannelLocations as a double


% --- Executes during object creation, after setting all properties.
function edit_ChannelLocations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ChannelLocations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function editImpedanceValues_Callback(hObject, eventdata, handles)
% hObject    handle to editImpedanceValues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editImpedanceValues as text
%        str2double(get(hObject,'String')) returns contents of editImpedanceValues as a double


% --- Executes during object creation, after setting all properties.
function editImpedanceValues_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editImpedanceValues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editElectrodeType_Callback(hObject, eventdata, handles)
% hObject    handle to editElectrodeType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editElectrodeType as text
%        str2double(get(hObject,'String')) returns contents of editElectrodeType as a double


% --- Executes during object creation, after setting all properties.
function editElectrodeType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editElectrodeType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editStimulusSource_Callback(hObject, eventdata, handles)
% hObject    handle to editStimulusSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editStimulusSource as text
%        str2double(get(hObject,'String')) returns contents of editStimulusSource as a double


% --- Executes during object creation, after setting all properties.
function editStimulusSource_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editStimulusSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMiscInfo_Callback(hObject, eventdata, handles)
% hObject    handle to editMiscInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMiscInfo as text
%        str2double(get(hObject,'String')) returns contents of editMiscInfo as a double


% --- Executes during object creation, after setting all properties.
function editMiscInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMiscInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
