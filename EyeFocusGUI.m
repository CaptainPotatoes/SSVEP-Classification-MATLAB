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

global countFar countMiddle countClose
% Choose default command line output for EyeFocusGUI
handles.output = hObject;
countFar = 1;
countMiddle = 1;
countClose = 1;

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
global myDevice Idx FarIdx MiddleIdx CloseIdx countFar countMiddle countClose
BioRadio_Name = 'EEG-SSVEP';
numEnabledBPChannels = double(myDevice.BioPotentialSignals.Count);

if numEnabledBPChannels == 0
    myDevice.Disconnect;
    BioRadioData = [];
    errordlg('No BioPotential Channels Programmed. Return to BioCapture to Configure.')
    return
end

sampleRate_BP = double(myDevice.BioPotentialSignals.SamplesPerSecond);
%Preallocating and setting up which area in the GUI the plot will go into
    % First two are raw data
    % Next two are data analysis features. 
numAxes = 4; 
axis_handles = zeros(1,numAxes);
for ch = 1:numAxes
    axis_handles(ch) = handles.(['axes',num2str(ch)]);
    if ch==1
        title([char(BioRadio_Name)]) 
    end
end

%Preallocating BPSignals

BioPotentialSignals = cell(1,numEnabledBPChannels);
Idx = cell(1,numEnabledBPChannels);
ButterFilt = cell(1,numEnabledBPChannels);
BufferFilt = cell (1,numEnabledBPChannels);
PSD = cell(1,numEnabledBPChannels);
%Filter
High=350*2/sampleRate_BP;
Low=450*2/sampleRate_BP; 
[b,a] = butter(3,[5*2/sampleRate_BP],'high');
%Notch Filter
NotchHigh = 126*2/sampleRate_BP;
NotchLow = 124*2/sampleRate_BP;
% [b1,a1] = butter(3,[NotchLow, NotchHigh],'stop');
if get(hObject,'Value') == 1
myDevice.StartAcquisition;
end
zold = 0;
plotWindow = 5;
plotGain_BP = 1;

while get(hObject,'Value') == 1
    pause(0.08)
    for ch = 1:numEnabledBPChannels

            BioPotentialSignals{ch} = [BioPotentialSignals{ch};myDevice.BioPotentialSignals.Item(ch-1).GetScaledValueArray.double'];
            Idx{ch} = 1:length(BioPotentialSignals{ch});
            %{
%             ButterFilt{ch}  = filtfilt(b,a,BioPotentialSignals{ch});
%             if length(BioPotentialSignals{1}) > 5000
%                 BufferFilt{1}  = buffer(BioPotentialSignals{ch},5000);
% 
%                 
% %                 for z = size(BufferFilt{1},2)
% %                    
% %                     if zold < z
% %                         if size(BufferFilt{1}(:,z),1) == 5000
% %                             zold = z;
% % 
% %     %                     Lo=size(BufferFilt{1},1);
% %     %                     NFFT = 2^nextpow2(Lo);
% %     %                     Yo = fft(BufferFilt{1}(:,z),NFFT)/Lo;
% %     %                     fo = sampleRate_BP/2*linspace(0,1,NFFT/2+1);
% % %                             [PSD{ch},f]=pwelch(BufferFilt{1}(:,z),[],[],[],sampleRate_BP);
% % % 
% % %                             plot(axis_handles(ch+2),f,pow2db(PSD{ch}))
% % %                             set(handles.(['axes',num2str(ch+2)]),'XLim',[0 f(end)])
% %                         end
% %                     end
% %                 end
%             end
            %}
            
            %Plot the Axes in the GUI
            if length(BioPotentialSignals{ch}) <= plotWindow*sampleRate_BP
                t = (0:(length(BioPotentialSignals{ch})-1))*(1/sampleRate_BP);
                plot(axis_handles(ch),t,plotGain_BP*BioPotentialSignals{ch})

    %             t2 = (0:(length(ButterFilt{ch})-1))*(1/sampleRate_BP);
    %             plot(axis_handles(ch+1),t2,plotGain_BP*ButterFilt{ch})

                set(handles.(['axes',num2str(ch)]),'XLim',[0 plotWindow]);
                set(get(handles.(['axes',num2str(ch)]), 'XLabel'), 'String', 'Time(s)')
                set(get(handles.(['axes',num2str(ch)]), 'YLabel'), 'String',  'mV')
                set(get(handles.(['axes',num2str(ch)]), 'Title'), 'String', 'EEG Bioradio')
            else
                if ch==1
                     t = ((length(BioPotentialSignals{ch})-(plotWindow*sampleRate_BP-1)):length(BioPotentialSignals{ch}))*(1/sampleRate_BP);
    %                  t2 = ((length(ButterFilt{ch})-(plotWindow*sampleRate_BP-1)):length(ButterFilt{ch}))*(1/sampleRate_BP);
                end
%             plot(axis_handles(ch+1),t2,plotGain_BP*ButterFilt{ch}(end-plotWindow*sampleRate_BP+1:end))
                plot(axis_handles(ch),t,plotGain_BP*BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end))
                set(handles.(['axes',num2str(ch)]),'XLim',[t(end)-plotWindow t(end)]);
                set(get(handles.(['axes',num2str(ch)]), 'XLabel'), 'String', 'Time(s)')
                set(get(handles.(['axes',num2str(ch)]), 'YLabel'), 'String',  'mV')
                set(get(handles.(['axes',num2str(ch)]), 'Title'), 'String', 'EEG Bioradio')
            end
            %% Todo: FFT and plot on axes #3
%             if length(BioPotentialSignals{ch}) > 500
%             plot(axis_handles(ch+1),lags{ch}(size(BufferFilt{1},2)-1,:)/sampleRate_BP,BPAutocorrelation{ch}(size(BufferFilt{1},2)-1,:))
%             end
            
   end
  
end

if get(hObject,'Value') == 0
    myDevice.StopAcquisition;
    BioRadioData = cell(1,4);
  
            BioRadioData{1,1} = BioPotentialSignals{1};
            BioRadioData{1,2} = BioPotentialSignals{2};
                %Analysis Stuff
            %BioRadioData{1,3} = ;
            %B
    
    assignin('base','Trial',BioRadioData)
%     countFar =1;
%     countMiddle =1;
%     countClose=1;
   
    
end

% assignin('base','FilteredSignal',ButterFilt)
% assignin('base','FFTData',PSD)
