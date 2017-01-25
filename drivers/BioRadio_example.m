%%
%
% load the BioRadio API using a MATLAB's .NET interface
%
current_dir = cd;

[ deviceManager , flag ] = load_API([current_dir '\BioRadioSDK.dll']);
% input = full path to api dll file
% outputs = deviceManager object, success flag
%
if ~flag % if API not successfully loaded, do not continue
    return
end
%
%
%%
%
% search for available sensors and select one
%
[ deviceName , macID , ok ] = BioRadio_Find( deviceManager );
% input = deviceManager object
% outputs = device name, macid, and flag if selection was canceled out
%

if ~ok %if no sensors selected, do not continue
    errordlg('Please select a BioRadio.')
    return
end

%
%
%%
%
% initialize BioRadio object
%
[ myDevice, flag ] = BioRadio_Connect ( deviceManager , macID , deviceName );
% input = deviceManager object, 64-bit mac address of BioRadio, and name of
% BioRadio
% outputs = BioRadio object, success flag for connection
%
if ~flag %if connection failed, do not continue
    return
end
% rather than search for BioRadios as outlined above, one can directly
% configure the code to initialize connection if the MAC ID is known
% e.g., myDevice = deviceManager.GetBluetoothDevice(int64(hex2dec('00a096388623')));
%
%
%%
%
% stream BioRadio data
%
BioRadioData = BioRadio_Stream2( myDevice , 10 , deviceName );
% inputs = BioRadio device handle, duration of data collection
% in seconds, string containing name of bioradio
% output = cell array of BioRadio data
BioRadioData1 = cell2mat(BioRadioData{1,1}(1:end,1));
BioRadioData2 = cell2mat(BioRadioData{1,1}(1:end,2));
BioRadioData3 = cell2mat(BioRadioData{1,1}(1:end,3));
%
%%
%
% disconnect from all sensors
%
BioRadio_Disconnect( myDevice )