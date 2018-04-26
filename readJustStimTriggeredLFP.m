function [LFPbySweep,Fs,completeSweeps]=readJustStimTriggeredLFP(dataDir,daqFileName,Fs,downSampleFactor,physChannels)
% Wrapper for daqread
% Reads the daq file containing the LFP and photodiode data in daqFileName
% 
% physChannel is the physiology channel of the rig (NI-DAQ channel number
% in Matlab)
% photoChannel is the photodiode channel of the rig (NI-DAQ channel number
% in Matlab)
% ledChannel is the LED channel of the rig (NI-DAQ channel number
% in Matlab)
%
% Fs is the sampling frequency of the data in daqFileName
%
% Returns an array of low-pass filtered (below lower_cutoff) and high-pass filtered
% (above upper_cutoff) data and the new sampling frequency of LFPbySweep
%
% Each row of LFPbySweep is the data from one sweep
% Different rows are different sweeps
% Different columns are different samples
%
% If bandPass==1: bandPass the data between lower_cutoff and upper_cutoff
% before returning; else return unfiltered LFP data
% 
% completeSweeps is a vector indicating which trials contain a complete
% acquisition sweep, that is, the entire sweep was acquired
% 
% Shorter sweeps (e.g., DaqController "hard-stopped" before sweep end) are
% padded with zeros
% Such shorter sweeps are indicated by zeros in completeSweeps; a value of
% 1 at index i in completeSweeps indicates a complete sweep in row # i
% of each returned array
%
% Choice: the rest of my analysis code either ignores incomplete trials or
% uses them (but note that averages will be affected by the trailing zeros)

%data=daqread(strcat(dataDir,'\',daqFileName),'Channels',[physChannel photoChannel ledChannel]);
data=daqread(strcat(dataDir,daqFileName),'Channels',physChannels);
% physData=data(:,physChannel);
% photoData=data(:,photoChannel);
% ledData=data(:,ledChannel);
for i=1:length(physChannels)
    physData{i}=data(:,i);
end

% Parse data into sweeps
for q=1:length(physChannels) 
    pData=physData{q};
    sweepSeps=find(isnan(pData));
    disp('Loading DAQ physiology data');
    if numel(sweepSeps)~=0
        sweepSepData=zeros(length(sweepSeps)+1,sweepSeps(1)-1);
        completeSweeps=zeros(length(sweepSeps)+1,1);
        sweepSepData(1,:)=pData(1:sweepSeps(1)-1)';
        completeSweeps(1)=1;
        for i=1:length(sweepSeps)
            if mod(i,100)==0
                disp('Loading DAQ physiology data');
            end
            if i==length(sweepSeps)
                if length(pData)-sweepSeps(end)==length(sweepSepData(1,:))
                    sweepSepData(i+1,:)=pData(sweepSeps(i)+1:end)';
                    completeSweeps(i+1)=1;
                else
                    sweepSepData(i+1,:)=[pData(sweepSeps(i)+1:end); zeros(length(sweepSepData(1,:))-(length(pData)-sweepSeps(end)),1)]';
                    completeSweeps(i+1)=0;
                end
            else
                sweepSepData(i+1,:)=pData(sweepSeps(i)+1:sweepSeps(i+1)-1)';
                completeSweeps(i+1)=1;
            end
        end
    else
        completeSweeps=1;
        sweepSepData=pData';
    end
    physData{q}=sweepSepData;
end

for q=1:length(physChannels)
    if downSampleFactor~=1
        currFs=Fs;
        [filtSweepSepData Fs]=downSamp(physData{q},currFs,downSampleFactor);
    else
        filtSweepSepData=physData{q};
    end
    physData{q}=filtSweepSepData;
end

LFPbySweep=physData;
