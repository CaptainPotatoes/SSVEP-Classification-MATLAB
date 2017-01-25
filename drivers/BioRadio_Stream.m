function  BioRadioData = BioRadio_Stream( myDevice , duration , BioRadio_Name )
% function  BioRadioData = BioRadio_Stream( myDevice , duration , BioRadio_Name )
% BioRadio_Stream streams data from the BioRadio and imports it into MATLAB.
%
% INPUTS:
% - myDevice is a handle to a BioRadio device object
% - duration is the data collection interval in seconds
% - BioRadio_name is string containing the BioRadio name
%
% OUTPUTS:
% - BioRadioData is a 3-element cell array of the data, where BioRadioData{1}
%   contains the BioPotentialSignals, BioRadioData{2} contains the
%   AuxiliarySignals, and BioRadioData{3} contains the PulseOxSignals. Each
%   of those cells contains a cell array where the number of cells
%   corresponds to the number of channels in that signal group. For
%   example, if 4 biopotentials are configured, BioRadioData{1} will be a 4
%   cell array, where each cell contains the data for a single channel.
%


numEnabledBPChannels = double(myDevice.BioPotentialSignals.Count);
numAuxChannels = double(myDevice.AuxiliarySignals.Count);
numPOxChannels = double(myDevice.PulseOxSignals.Count);

if numEnabledBPChannels == 0
    myDevice.Disconnect;
    BioRadioData = [];
    errordlg('No BioPotential Channels Programmed. Return to BioCapture to Configure.')
    return
end

sampleRate_BP = double(myDevice.BioPotentialSignals.SamplesPerSecond);
sampleRate_Pod = 250;

figure
axis_handles = zeros(1,numEnabledBPChannels+numAuxChannels+numPOxChannels);
for ch = 1:numEnabledBPChannels
    axis_handles(ch) = subplot(length(axis_handles),1,ch);
    if ch==1
        title([char(BioRadio_Name)])
    end
    ylabel([char(myDevice.BioPotentialSignals.Item(ch-1).Name) ' (V)']);
%     ylim(10*double([myDevice.BioPotentialSignals.Item(ch-1).MinValue myDevice.BioPotentialSignals.Item(ch-1).MaxValue]))
    hold on
end
for ch = 1:numAuxChannels
    axis_handles(ch+numEnabledBPChannels) = subplot(length(axis_handles),1,ch+numEnabledBPChannels);
    ylabel(char(myDevice.AuxiliarySignals.Item(ch-1).Name));
    ylim(double([myDevice.AuxiliarySignals.Item(ch-1).MinValue myDevice.AuxiliarySignals.Item(ch-1).MaxValue]))
    hold on
end
for ch = 1:numPOxChannels
    axis_handles(ch+numEnabledBPChannels+numAuxChannels) = subplot(length(axis_handles),1,ch+numEnabledBPChannels+numAuxChannels);
    ylabel(char(myDevice.PulseOxSignals.Item(ch-1).Name));
    ylim(double([myDevice.PulseOxSignals.Item(ch-1).MinValue myDevice.PulseOxSignals.Item(ch-1).MaxValue]))
    hold on
end
xlabel('Time (s)')

linkaxes(axis_handles,'x')



BioPotentialSignals = cell(1,numEnabledBPChannels);

AuxiliarySignals = cell(1,numAuxChannels);

PulseOxSignals = cell(1,numPOxChannels);

myDevice.StartAcquisition;

plotWindow = 5;

plotGain_BP = 1;

elapsedTime = 0;
tic;

while elapsedTime < duration
    pause(0.08)
    for ch = 1:numEnabledBPChannels
        BioPotentialSignals{ch} = [BioPotentialSignals{ch};myDevice.BioPotentialSignals.Item(ch-1).GetScaledValueArray.double'];
        if length(BioPotentialSignals{ch}) <= plotWindow*sampleRate_BP
            cla(axis_handles(ch))
            t = (0:(length(BioPotentialSignals{ch})-1))*(1/sampleRate_BP);
            plot(axis_handles(ch),t,plotGain_BP*BioPotentialSignals{ch});
            xlim([0 plotWindow])
        else
            if ch==1
                t = ((length(BioPotentialSignals{ch})-(plotWindow*sampleRate_BP-1)):length(BioPotentialSignals{ch}))*(1/sampleRate_BP);
            end
            cla(axis_handles(ch))
            plot(axis_handles(ch),t,plotGain_BP*BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end));
            xlim([t(end)-plotWindow t(end)])
        end
    end
    
    for ch = 1:numAuxChannels
        AuxiliarySignals{ch} = [AuxiliarySignals{ch};myDevice.AuxiliarySignals.Item(ch-1).GetScaledValueArray.double'];
        if length(AuxiliarySignals{ch}) <= plotWindow*sampleRate_Pod
            cla(axis_handles(ch+numEnabledBPChannels))
            t = (0:(length(AuxiliarySignals{ch})-1))*(1/sampleRate_Pod);
            plot(axis_handles(ch+numEnabledBPChannels),t,AuxiliarySignals{ch});
            xlim([0 plotWindow])
        else
            if ch==1
                t_pod = ((length(AuxiliarySignals{ch})-(plotWindow*sampleRate_Pod-1)):length(AuxiliarySignals{ch}))*(1/sampleRate_Pod);
            end
            cla(axis_handles(ch+numEnabledBPChannels))
            plot(axis_handles(ch+numEnabledBPChannels),t_pod,AuxiliarySignals{ch}(end-plotWindow*sampleRate_Pod+1:end));
            xlim([t(end)-plotWindow t(end)])
        end
        hold on
    end
    
    for ch = 1:numPOxChannels
        PulseOxSignals{ch} = [PulseOxSignals{ch};myDevice.PulseOxSignals.Item(ch-1).GetScaledValueArray.double'];
        if length(PulseOxSignals{ch}) <= plotWindow*sampleRate_Pod
            cla(axis_handles(ch+numEnabledBPChannels+numAuxChannels))
            t = (0:(length(PulseOxSignals{ch})-1))*(1/sampleRate_Pod);
            plot(axis_handles(ch+numEnabledBPChannels+numAuxChannels),t,PulseOxSignals{ch});
            xlim([0 plotWindow])
        else
            if ch==1
                t_pod = ((length(PulseOxSignals{ch})-(plotWindow*sampleRate_Pod-1)):length(PulseOxSignals{ch}))*(1/sampleRate_Pod);
            end
            cla(axis_handles(ch+numEnabledBPChannels+numAuxChannels))
            plot(axis_handles(ch+numEnabledBPChannels+numAuxChannels),t_pod,PulseOxSignals{ch}(end-plotWindow*sampleRate_Pod+1:end));
            xlim([t(end)-plotWindow t(end)])
        end
    end
    
    elapsedTime = elapsedTime + toc;
    tic;
end

myDevice.StopAcquisition;

for ch = 1:numEnabledBPChannels
    BioPotentialSignals{ch} = [BioPotentialSignals{ch};myDevice.BioPotentialSignals.Item(ch-1).GetScaledValueArray.double'];
    t = ((length(BioPotentialSignals{ch})-(plotWindow*sampleRate_BP-1)):length(BioPotentialSignals{ch}))*(1/sampleRate_BP);
    cla(axis_handles(ch))
    plot(axis_handles(ch),t,plotGain_BP*BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end));
    xlim([t(end)-plotWindow t(end)])
end

for ch = 1:numAuxChannels
    AuxiliarySignals{ch} = [AuxiliarySignals{ch};myDevice.AuxiliarySignals.Item(ch-1).GetScaledValueArray.double'];
    t = ((length(AuxiliarySignals{ch})-(plotWindow*sampleRate_Pod-1)):length(AuxiliarySignals{ch}))*(1/sampleRate_Pod);
    cla(axis_handles(ch+numEnabledBPChannels))
    plot(axis_handles(ch+numEnabledBPChannels),t,AuxiliarySignals{ch}(end-plotWindow*sampleRate_Pod+1:end));
    xlim([t(end)-plotWindow t(end)])
end


for ch = 1:numPOxChannels
    PulseOxSignals{ch} = [PulseOxSignals{ch};myDevice.PulseOxSignals.Item(ch-1).GetScaledValueArray.double'];
    t = ((length(PulseOxSignals{ch})-(plotWindow*sampleRate_Pod-1)):length(PulseOxSignals{ch}))*(1/sampleRate_Pod);
    cla(axis_handles(ch+numEnabledBPChannels+numAuxChannels))
    plot(axis_handles(ch+numEnabledBPChannels+numAuxChannels),t,PulseOxSignals{ch}(end-plotWindow*sampleRate_Pod+1:end));
    xlim([t(end)-plotWindow t(end)])
end

BioRadioData = cell(1,3);
BioRadioData{1} = BioPotentialSignals;
BioRadioData{2} = AuxiliarySignals;
BioRadioData{3} = PulseOxSignals;

end