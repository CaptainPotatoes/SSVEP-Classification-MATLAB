function BioRadio_Disconnect ( myDevice )
% function BioRadio_Disconnect ( myDevice )
% BioRadio_Disconnect terminates the connection with the BioRadio
% 
% INPUTS:
% - myDevice is a handle to a BioRadio device object

myDevice.Disconnect;