function varargout = WheelchairManualControls(varargin)
% WHEELCHAIRMANUALCONTROLS MATLAB code for WheelchairManualControls.fig
%      WHEELCHAIRMANUALCONTROLS, by itself, creates a new WHEELCHAIRMANUALCONTROLS or raises the existing
%      singleton*.inde
%
%      H = WHEELCHAIRMANUALCONTROLS returns the handle to a new WHEELCHAIRMANUALCONTROLS or the handle to
%      the existing singleton*.
%
%      WHEELCHAIRMANUALCONTROLS('CALLBACK',hObject, eventData,handles,...) calls the local
%      function named CALLBACK in WHEELCHAIRMANUALCONTROLS.M with the given input arguments.
%
%      WHEELCHAIRMANUALCONTROLS('Property','Value',...) creates a new WHEELCHAIRMANUALCONTROLS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WheelchairManualControls_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WheelchairManualControls_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WheelchairManualControls

% Last Modified by GUIDE v2.5 01-Jun-2016 13:42:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WheelchairManualControls_OpeningFcn, ...
                   'gui_OutputFcn',  @WheelchairManualControls_OutputFcn, ...
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


% --- Executes just before WheelchairManualControls is made visible.
function WheelchairManualControls_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WheelchairManualControls (see VARARGIN)

% Choose default command line output for WheelchairManualControls
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes WheelchairManualControls wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WheelchairManualControls_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in UpButton.
function UpButton_Callback(hObject, eventdata, handles)
global b1
fwrite(b1,'w')
% hObject    handle to UpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in LeftButton.
function LeftButton_Callback(hObject, eventdata, handles)
global b1
fwrite(b1,'d')
% hObject    handle to LeftButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in DownButton.
function DownButton_Callback(hObject, eventdata, handles)
global b1
fwrite(b1,'s')

% hObject    handle to DownButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in RightButton.
function RightButton_Callback(hObject, eventdata, handles)
global b1
fwrite(b1,'a')

% hObject    handle to RightButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btButton.
function btButton_Callback(hObject, eventdata, handles)
% hObject    handle to btButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global b1
BluetoothButton = get(hObject,'Value');
if BluetoothButton == 1
          disp(BluetoothButton)
b1 = Bluetooth('Adafruit EZ-Link 9580',1);
%pause(3)
fopen(b1)
elseif BluetoothButton == 0
         disp(BluetoothButton)
     b1 = Bluetooth('Adafruit EZ-Link 9580',1);

fclose(b1)

end

% Hint: get(hObject,'Value') returns toggle state of btButton

% --- Executes on button press in safetyButton.
function safetyButton_Callback(hObject, eventdata, handles)
% hObject    handle to safetyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SafetyButtonToggleState = get(hObject,'Value');

if SafetyButtonToggleState == 1
        disp(SafetyButtonToggleState)

elseif SafetyButtonToggleState == 0
        disp(SafetyButtonToggleState)

end


% Hint: get(hObject,'Value') returns toggle state of safetyButton
