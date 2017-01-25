function [ deviceManager , flag ] = load_API( filepath )
% function [ deviceManager , flag ] = load_API( filepath )
% load_API.m makes the BioRadio .NET assembly visible to MATLAB
% 
% INPUTS:
% - filepath is the file path to the API dll. If it is not provided as
%   an input, the user is prompted to select a file. 
% 
% OUTPUTS:
% - flag that is a boolean that is true if the assembly successfully 
%   loaded and false otherwise.
% - deviceManager is a handle to BioRadio device manager object
%

switch nargin
    case 0 % if filepath not provided prompt user to select dll
        [filename,pathname]=uigetfile('*.dll','Select the API dll file.');
        try
            asmInfo = NET.addAssembly([pathname filename]); %load API
            deviceManager = GLNeuroTech.Devices.BioRadio.BioRadioDeviceManager; %create device manager object
            flag = true;
        catch
            errordlg('Invalid file selection. Please try again')
            flag = false;
            deviceManager = [];
            return
        end
    case 1 %if filepath provided, use user input to select dll
        try
            asmInfo = NET.addAssembly(filepath); %load API
            deviceManager = GLNeuroTech.Devices.BioRadio.BioRadioDeviceManager; %create device manager object
            flag = true;
            return
        catch
            errordlg('Invalid file selection. Please try again')
            flag = false;
            deviceManager = [];
            return
        end
    otherwise % if more than two inputs provided, error out
        errordlg('Invalid file selection. Please try again')
        flag = false;
        deviceManager = [];
        return
end


end

