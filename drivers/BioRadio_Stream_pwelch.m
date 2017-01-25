function  BioRadioData = BioRadio_Stream2( myDevice , duration , BioRadio_Name )
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

if numEnabledBPChannels == 0
    myDevice.Disconnect;
    BioRadioData = [];
    errordlg('No BioPotential Channels Programmed. Return to BioCapture to Configure.')
    return
end

sampleRate_BP = double(myDevice.BioPotentialSignals.SamplesPerSecond);


figure
axis_handles = zeros(1,numEnabledBPChannels);
for ch = 1:numEnabledBPChannels
    axis_handles(ch) = subplot(length(axis_handles),1,ch);
    if ch==1
        title([char(BioRadio_Name)])
    end
    ylabel([char(myDevice.BioPotentialSignals.Item(ch-1).Name) ' (V)']);
%     ylim(10*double([myDevice.BioPotentialSignals.Item(ch-1).MinValue myDevice.BioPotentialSignals.Item(ch-1).MaxValue]))
    hold on
end
xlabel('Time (s)')
linkaxes(axis_handles,'x')

figure
axis_handles2 = zeros(1,numEnabledBPChannels);
for ch = 1:numEnabledBPChannels
    axis_handles(ch) = subplot(length(axis_handles),1,ch);
    if ch==1
        title([char(BioRadio_Name)])
    end
    ylabel([char(myDevice.BioPotentialSignals.Item(ch-1).Name) 'Mag']);
%     ylim(10*double([myDevice.BioPotentialSignals.Item(ch-1).MinValue myDevice.BioPotentialSignals.Item(ch-1).MaxValue]))
    hold on
end
xlabel('Frequency (hz)')





BioPotentialSignals = cell(1,numEnabledBPChannels);
p = cell(1,numEnabledBPChannels);
f= cell(1,numEnabledBPChannels);
% BioPotentialSignalsFiltered = cell(2,numEnabledBPChannels);
% rads=.1*2/1000;
% pass=4*2/1000;
% [b,a] = butter(2,[ pass ,rads],'bandpass');


myDevice.StartAcquisition;

plotWindow = 5;

plotGain_BP = 1;

elapsedTime = 0;
tic;

while elapsedTime < duration
    pause(0.08)
    for ch = 1:numEnabledBPChannels

        BioPotentialSignals{ch} = [BioPotentialSignals{ch};myDevice.BioPotentialSignals.Item(ch-1).GetScaledValueArray.double'];
%         BioPotentialSignalsFiltered{ch} = [filtered(BioPotentialSignals{ch};myDevice.BioPotentialSignals.Item(ch-1).GetScaledValueArray.double'];
      
   
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
  
   for ch = 1:numEnabledBPChannels
      [p{ch},f{ch}]=pwelch(BioPotentialSignals{ch},1000,1000,sampleRate_BP);
      if length(BioPotentialSignals{ch})<=plotWindow*sampleRate_BP
          cla(axis_handles2(ch))
          plot(axis_handles2(ch),f{ch},p{ch})
      else
          cla(axis_handles2(ch))
          plot(axis_handles2(ch),f{ch},p{ch}(end-plotWindow*sampleRate_BP+1:end));
      end
   end   
       
    elapsedTime = elapsedTime + toc;
    tic;
end

myDevice.StopAcquisition;

for ch = 1:numEnabledBPChannels
%     BioPotentialSignals{ch} = filter(b,a,BioPotentialSignals{ch})+ BioPotentialSignals{ch};
    BioPotentialSignals{ch} = [BioPotentialSignals{ch};myDevice.BioPotentialSignals.Item(ch-1).GetScaledValueArray.double'];
    t = ((length(BioPotentialSignals{ch})-(plotWindow*sampleRate_BP-1)):length(BioPotentialSignals{ch}))*(1/sampleRate_BP);
    cla(axis_handles(ch))
    plot(axis_handles(ch),t,plotGain_BP*BioPotentialSignals{ch}(end-plotWindow*sampleRate_BP+1:end));
    xlim([t(end)-plotWindow t(end)])
end



BioRadioData = cell(1,1);
BioRadioData{1}= BioPotentialSignals;
BioRadioData{2}= BioPotentialSignalsfiltered;
assignin('base', 'BioRadioData{2}', BioRadioData{2})


end


