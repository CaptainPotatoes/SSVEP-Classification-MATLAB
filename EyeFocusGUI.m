function varargout = EyeFocusGUI(varargin)
% EYEFOCUSGUI MATLAB code for EyeFocusGUI.fig
%      EYEFOCUSGUI, by itself, creates a new EYEFOCUSGUI or raises the existing
%      singleton*.
%
%      H = EYEFOCUSGUI returns the handle to a new EYEFOCUSGUI or the handle to
%      the existing singleton*.
%
%      EYEFOCUSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EYEFOCUSGUI.M with the given input arguments.
%
%      EYEFOCUSGUI('Property','Value',...) creates a new EYEFOCUSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EyeFocusGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EyeFocusGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help EyeFocusGUI
% Last Modified by GUIDE v2.5 29-Jan-2017 12:44:30
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EyeFocusGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @EyeFocusGUI_OutputFcn, ...
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

% --- Executes just before EyeFocusGUI is made visible.
function EyeFocusGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EyeFocusGUI (see VARARGIN)

% global countFar countMiddle countClose
% Choose default command line output for EyeFocusGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EyeFocusGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = EyeFocusGUI_OutputFcn(hObject, eventdata, handles) 
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
count = 0;

if get(hObject,'Value') && count == 0

% load the BioRadio API using a MATLAB's .NET interface
[ deviceManager , flag ] = load_API(['C:\Users\mahmoodms\Dropbox\Public\_VCU\Yeo Lab\_SSVEP\_MATLAB-SSVEP-Classification\BioRadioSDK.dll']);
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
BioRadio_Name = 'EEG-SSVEP';
numEnabledBPChannels = double(myDevice.BioPotentialSignals.Count);

if numEnabledBPChannels == 0
    myDevice.Disconnect;
    RawBioRadioData = [];
    errordlg('No BioPotential Channels Programmed. Return to BioCapture to Configure.')
    return
end

sampleRate_BP = double(myDevice.BioPotentialSignals.SamplesPerSecond);
%Preallocating and setting up which area in the GUI the plot will go into
    % First two are raw data
    % Next two are data analysis features. 
numAxes = 8; 
axis_handles = zeros(1,numAxes);
for ch = 1:numAxes
    axis_handles(ch) = handles.(['axes',num2str(ch)]);
    %{
    if ch==1
        title([char(BioRadio_Name)]) 
    end
    %}
end

%Preallocating BPSignals

BioPotentialSignals = cell(1,numEnabledBPChannels);
Idx = cell(1,numEnabledBPChannels);
if get(hObject,'Value') == 1
myDevice.StartAcquisition;
end
plotWindow = 5;
plotGain_BP = 1;
% fft_len = plotWindow*sampleRate_BP;
fft_len = 2^(nextpow2(plotWindow*sampleRate_BP));
% dBmax = 100;
spect_1 = handles.axes5;
spect_2 = handles.axes6;
max_freq = 50;
while get(hObject,'Value') == 1
    pause(0.08)
    for ch = 1:numEnabledBPChannels
            BioPotentialSignals{ch} = [BioPotentialSignals{ch};myDevice.BioPotentialSignals.Item(ch-1).GetScaledValueArray.double'];
            Idx{ch} = 1:length(BioPotentialSignals{ch});            
            %Plot the Axes in the GUI
            if length(BioPotentialSignals{ch}) <= plotWindow*sampleRate_BP
                t = (0:(length(BioPotentialSignals{ch})-1))*(1/sampleRate_BP);
                plot(axis_handles(ch),t,plotGain_BP*BioPotentialSignals{ch})
                set(handles.(['axes',num2str(ch)]),'XLim',[0 plotWindow]);
                set(get(handles.(['axes',num2str(ch)]), 'XLabel'), 'String', 'Time(s)')
                set(get(handles.(['axes',num2str(ch)]), 'YLabel'), 'String',  'mV')
                if ch==1
                    set(get(handles.(['axes',num2str(ch)]), 'Title'), 'String', 'Fp1')
                    if length(BioPotentialSignals{ch})>sampleRate_BP
                        fft_len_short = 2^(nextpow2(length(BioPotentialSignals{ch})));
                        %FILTER:
                        fp1_fft = fft(eeg_h_fcn3_50(BioPotentialSignals{ch},sampleRate_BP));
                        P2 = abs(fp1_fft/fft_len_short);
                        P1 = P2(1:fft_len_short/2+1);
                        P1(2:end-1) = 2*P1(2:end-1);
                        f = sampleRate_BP*(0:(fft_len_short/2))/fft_len_short;
                        plot(axis_handles(3),f,P1);
                        set(handles.(['axes',num2str(3)]),'XLim',[1 100]);
                        set(get(handles.(['axes',num2str(3)]), 'XLabel'), 'String', 'f (Hz)')
                        set(get(handles.(['axes',num2str(3)]), 'YLabel'), 'String', '|P1(f)|')
                        set(get(handles.(['axes',num2str(3)]), 'Title'), 'String', 'FFT(Fp1)')
                    end
                elseif ch==2
                    set(get(handles.(['axes',num2str(ch)]),'Title'),'String','Fp2')
                    if length(BioPotentialSignals{ch})>sampleRate_BP
                        fft_len_short = 2^(nextpow2(length(BioPotentialSignals{ch})));
                        %FILTER:
                        fp2_fft = fft(eeg_h_fcn3_50(BioPotentialSignals{ch},sampleRate_BP));
                        P2 = abs(fp2_fft/fft_len_short);
                        P1 = P2(1:fft_len_short/2+1);
                        P1(2:end-1) = 2*P1(2:end-1);
                        f = sampleRate_BP*(0:(fft_len_short/2))/fft_len_short;
                        plot(axis_handles(4),f,P1);
                        set(handles.(['axes',num2str(4)]),'XLim',[1 100]);
                        set(get(handles.(['axes',num2str(4)]), 'XLabel'), 'String', 'f (Hz)')
                        set(get(handles.(['axes',num2str(4)]), 'YLabel'), 'String', '|P1(f)|')
                        set(get(handles.(['axes',num2str(4)]), 'Title'), 'String', 'FFT(Fp2)')
                    end
                end
            else %once plot window is exceeded:
                if ch==1
                     t = ((length(BioPotentialSignals{ch})-(plotWindow*sampleRate_BP-1)):length(BioPotentialSignals{ch}))*(1/sampleRate_BP);
                end
                plot(axis_handles(ch),t,plotGain_BP*BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end))
                set(handles.(['axes',num2str(ch)]),'XLim',[t(end)-plotWindow t(end)]);
                set(get(handles.(['axes',num2str(ch)]), 'XLabel'), 'String', 'Time(s)')
                set(get(handles.(['axes',num2str(ch)]), 'YLabel'), 'String',  'mV')
                if ch==1
%                     set(get(handles.(['axes',num2str(ch)]), 'Title'), 'String', 'Fp1')
                    % FFT:
                    fp1_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end);
                    fp1_data_filtered = eeg_h_fcn3_50(fp1_data_unfilt, sampleRate_BP);
                    fp1_fft = fft(fp1_data_filtered);
                    P2 = abs(fp1_fft/fft_len);
                    P1 = P2(1:fft_len/2+1);
                    P1(2:end-1) = 2*P1(2:end-1);
                    f = sampleRate_BP*(0:(fft_len/2))/fft_len;
                    plot(axis_handles(3),f,P1);
                    set(handles.(['axes',num2str(3)]),'XLim',[1 100]);
                    set(get(handles.(['axes',num2str(3)]), 'XLabel'), 'String', 'f (Hz)')
                    set(get(handles.(['axes',num2str(3)]), 'YLabel'), 'String', '|P1(f)|')
                    set(get(handles.(['axes',num2str(3)]), 'Title'), 'String', 'FFT(Fp1)')
                    % Spect:
                    [S, Fspect, T, P] = spectrogram(fp1_data_filtered, 5*sampleRate_BP,4*sampleRate_BP,10*sampleRate_BP,sampleRate_BP);
                    h1 = imagesc(spect_1, T, Fspect(Fspect<100 & Fspect>0), 10*log10(P(Fspect<100 & Fspect>0,:)));
                    set(spect_1,'YDir','normal')
                    cb = colorbar(spect_1);
                    ylabel(cb, 'Power (db)')
                    colormap(spect_1,jet)
                    set(get(handles.(['axes',num2str(5)]), 'XLabel'), 'String', 'Time (s)')
                    set(get(handles.(['axes',num2str(5)]), 'YLabel'), 'String', 'Frequency (Hz)')
                    set(get(handles.(['axes',num2str(5)]), 'Title'), 'String', 'Spectrogram (Fp1)')
                    %Pwelch
                    [Pxx, F] = pwelch(fp1_data_filtered,[],[],250);
                    plot(axis_handles(7), Pxx); 
                    set(handles.(['axes',num2str(7)]),'XLim',[0 100]);
                    set(get(handles.(['axes',num2str(7)]), 'XLabel'), 'String', 'Frequency (Hz)')
                    set(get(handles.(['axes',num2str(7)]), 'YLabel'), 'String', 'Power (dB)')
                    set(get(handles.(['axes',num2str(7)]), 'Title'), 'String', 'Pwelch (Fp1)')
                elseif ch==2
                    fp2_data_unfilt = BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end);
                    fp2_data_filtered = eeg_h_fcn3_50(fp2_data_unfilt, sampleRate_BP);
                    fp2_fft = fft(fp2_data_filtered);
                    P2 = abs(fp2_fft/fft_len);
                    P1 = P2(1:fft_len/2+1);
                    P1(2:end-1) = 2*P1(2:end-1);
                    f = sampleRate_BP*(0:(fft_len/2))/fft_len;
                    plot(axis_handles(4),f,P1);
                    set(handles.(['axes',num2str(4)]),'XLim',[1 100]);
                    set(get(handles.(['axes',num2str(4)]), 'XLabel'), 'String', 'f (Hz)')
                    set(get(handles.(['axes',num2str(4)]), 'YLabel'), 'String', '|P1(f)|')
                    set(get(handles.(['axes',num2str(4)]), 'Title'), 'String', 'FFT(Fp2)')
                    % Spect:
                    [S, Fspect, T, P] = spectrogram(fp2_data_filtered, 5*sampleRate_BP,4*sampleRate_BP,10*sampleRate_BP,sampleRate_BP);
                    h2 = imagesc(spect_2, T, Fspect(Fspect<100 & Fspect>0), 10*log10(P(Fspect<100 & Fspect>0,:)));
                    set(spect_2,'YDir','normal') %gca
                    cb2 = colorbar(spect_2);
                    ylabel(cb2, 'Power (db)')
                    colormap(spect_2,jet)
                    set(get(handles.(['axes',num2str(6)]), 'XLabel'), 'String', 'Time (s)')
                    set(get(handles.(['axes',num2str(6)]), 'YLabel'), 'String', 'Frequency (Hz)')
                    set(get(handles.(['axes',num2str(6)]), 'Title'), 'String', 'Spectrogram (Fp2)')
                    %Pwelch
                    [Pxx, F] = pwelch(fp2_data_filtered,[],[],250);
                    plot(axis_handles(8), Pxx); 
                    set(handles.(['axes',num2str(8)]),'XLim',[0 100]);
                    set(get(handles.(['axes',num2str(8)]), 'XLabel'), 'String', 'Frequency (Hz)')
                    set(get(handles.(['axes',num2str(8)]), 'YLabel'), 'String', 'Power (dB)')
                    set(get(handles.(['axes',num2str(8)]), 'Title'), 'String', 'Pwelch (Fp1)')
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
    RawBioRadioData = cell(1,2);
            RawBioRadioData{1,1} = BioPotentialSignals{1};
            RawBioRadioData{1,2} = BioPotentialSignals{2};
    assignin('base','Trial',RawBioRadioData)
    %% Todo: Change into button function.
    c=clock;
    filename = ['RawBioRadioData_',num2str(c(2)),'-',num2str(c(3)),'-',num2str(c(1)),'_',num2str(c(4)),'.',num2str(c(5)),',',num2str(c(6)),'.xlsx'];
    l1 = length(BioPotentialSignals{1})
    l2 = length(BioPotentialSignals{2})
    if l2~=l1 %if not the same size
        if l2>l1
            xlswrite(filename, [BioPotentialSignals{1},BioPotentialSignals{2}(1:l1)]);
        else %l1>l2
            xlswrite(filename, [BioPotentialSignals{1}(1:l2),BioPotentialSignals{2}]);
        end
    else
        xlswrite(filename, [BioPotentialSignals{1},BioPotentialSignals{2}]);
    end
end

% assignin('base','FilteredSignal',ButterFilt)
% assignin('base','FFTData',PSD)
