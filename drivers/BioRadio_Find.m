function [BioRadioName, macID, ok] = BioRadio_Find( deviceManager )
% function [name, macID, ok] = BioRadio_Find( deviceManager )
% BioRadio_Find searches for BioRadios in the vicinity and prompts the
% user to select one for connection
%
% INPUTS:
% - deviceManager is a device manager object created using the BioRadio API
%
% OUTPUTS:
% - name is a string containing the name of the selected BioRadio
% - macID the 64-bit mac address of the selected BioRadio
% - ok is a boolean flag indicating if user cancelled out of the selection 
%   menu and did not pick a BioRadio
%

dlghandle = helpdlg('Please wait... search in progress.','Searching for Available BioRadios');

try
    BioRadioList = deviceManager.DiscoverBluetoothDevices; %search for BioRadios
    numavail = BioRadioList.Length; 
catch
    numavail = 0;
end

close(dlghandle)

if numavail<1
    BioRadioName = [];
    macID = [];
    ok = [];
    return
end

availableBioRadios = cell(numavail,1);
macIDs = cell(numavail,1);

for i=1:numavail
    availableBioRadios{i} = char(BioRadioList(i).DeviceId); % pull BioRadio name from list
    macIDs{i} = hex2dec(char(BioRadioList(i).MacId));% pull corresponding mac id from list
end

[selection, ok] = listdlg('PromptString','Select a BioRadio:',...
    'SelectionMode','single','ListString',availableBioRadios); % prompt user to select a BioRadio

if ok==0
    BioRadioName = [];
    macID = [];
    ok = [];
    return
else
    BioRadioName = availableBioRadios{selection};
    macID = macIDs{selection};
end

end