function [ myDevice, flag ] = BioRadio_Connect ( deviceManager , macID , BioRadioName )
% function [ myDevice, flag ] = BioRadio_Connect ( deviceManager , macID , BioRadioName )
% BioRadio_Connect initializes the BioRadio objects and connects
% 
% INPUTS:
% - deviceManager is a handle to a motion sensor device manager object
% - macID is a 64-bit MAC ID (e.g., macID = int64(hex2dec(‘00A096388623’)))
% - BioRadioName is a string containing the BioRadio name. These are 2 character identifiers,
%   but alternatively, could be descriptive labels

%
% OUTPUTS:
% - myDevice is a handle to a BioRadio device object
% - flag is a logical success flag that is false if connection fails
%


myDevice = []; % create holder for motion sensor objects
flag = false;

try
    myDevice = deviceManager.GetBluetoothDevice(macID); % instantiate motion sensor object
catch
    errordlg(['Failed to connect to ' BioRadioName '.'])
    return
end

flag = true; % successfully connected to all sensors

end