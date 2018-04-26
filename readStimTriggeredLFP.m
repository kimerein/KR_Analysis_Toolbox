function [LFPbySweep,Fs,photoDatabySweep,ledDatabySweep,completeSweeps]=readStimTriggeredLFP(daqFileName,Fs,bandPass,lower_cutoff,upper_cutoff,downSampleFactor,physChannel,photoChannel,ledChannel)
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

global dataDir

%data=daqread(strcat(dataDir,'\',daqFileName),'Channels',[physChannel photoChannel ledChannel]);
data=daqread(strcat(dataDir,daqFileName),'Channels',[physChannel photoChannel ledChannel]);
% physData=data(:,physChannel);
% photoData=data(:,photoChannel);
% ledData=data(:,ledChannel);
physData=data(:,1);
photoData=data(:,2);
ledData=data(:,3);


% Parse data into sweeps
sweepSeps=find(isnan(physData));
disp('Loading DAQ physiology data');
if numel(sweepSeps)~=0
    sweepSepData=zeros(length(sweepSeps)+1,sweepSeps(1)-1);
    completeSweeps=zeros(length(sweepSeps)+1,1);
    sweepSepData(1,:)=physData(1:sweepSeps(1)-1)';
    completeSweeps(1)=1;
    for i=1:length(sweepSeps)
        if mod(i,100)==0
            disp('Loading DAQ physiology data');
        end
        if i==length(sweepSeps)
            if length(physData)-sweepSeps(end)==length(sweepSepData(1,:))
                sweepSepData(i+1,:)=physData(sweepSeps(i)+1:end)';
                completeSweeps(i+1)=1;
            else
                sweepSepData(i+1,:)=[physData(sweepSeps(i)+1:end); zeros(length(sweepSepData(1,:))-(length(physData)-sweepSeps(end)),1)]';
                completeSweeps(i+1)=0;
            end
        else
            sweepSepData(i+1,:)=physData(sweepSeps(i)+1:sweepSeps(i+1)-1)';
            completeSweeps(i+1)=1;
        end
    end
else
    completeSweeps=1;
    sweepSepData=physData';
end

% Parse photoData into sweeps
sweepSeps=find(isnan(photoData));
disp('Loading DAQ photodiode data');
if numel(sweepSeps)~=0
    photoSepData=zeros(length(sweepSeps)+1,sweepSeps(1)-1);
    photoSepData(1,:)=photoData(1:sweepSeps(1)-1)';
    for i=1:length(sweepSeps)
        if mod(i,100)==0
            disp('Loading DAQ photodiode data');
        end
        if i==length(sweepSeps)
            if length(photoData)-sweepSeps(end)==length(photoSepData(1,:))
                photoSepData(i+1,:)=photoData(sweepSeps(i)+1:end)';
            else
                photoSepData(i+1,:)=[photoData(sweepSeps(i)+1:end); zeros(length(photoSepData(1,:))-(length(photoData)-sweepSeps(end)),1)]';
            end
        else
            photoSepData(i+1,:)=photoData(sweepSeps(i)+1:sweepSeps(i+1)-1)';
        end
    end
else
    photoSepData=photoData';
end

% Parse ledData into sweeps
sweepSeps=find(isnan(ledData));
disp('Loading DAQ LED data');
if numel(sweepSeps)~=0
    ledSepData=zeros(length(sweepSeps)+1,sweepSeps(1)-1);
    ledSepData(1,:)=ledData(1:sweepSeps(1)-1)';
    for i=1:length(sweepSeps)
        if mod(i,100)==0
            disp('Loading DAQ LED data');
        end
        if i==length(sweepSeps)
            if length(ledData)-sweepSeps(end)==length(ledSepData(1,:))
                ledSepData(i+1,:)=ledData(sweepSeps(i)+1:end)';
            else
                ledSepData(i+1,:)=[ledData(sweepSeps(i)+1:end); zeros(length(ledSepData(1,:))-(length(ledData)-sweepSeps(end)),1)]';
            end
        else
            ledSepData(i+1,:)=ledData(sweepSeps(i)+1:sweepSeps(i+1)-1)';
        end
    end
else
    ledSepData=ledData';
end

if downSampleFactor~=1
    currFs=Fs;
    [filtSweepSepData Fs]=downSamp(sweepSepData,currFs,downSampleFactor);
    [photoDatabySweep pFs]=downSamp(photoSepData,currFs,downSampleFactor);
    [ledDatabySweep lFs]=downSamp(ledSepData,currFs,downSampleFactor);
else
    filtSweepSepData=sweepSepData;
    photoDatabySweep=photoSepData;
    ledDatabySweep=ledSepData;
end

if bandPass==1
    % Low-pass-filter data
    disp('LP-filtering');
    filtSweepSepData=fftFilter(filtSweepSepData',Fs,upper_cutoff,1);
    filtSweepSepData=filtSweepSepData';
    disp('Done LP-filtering');
    
    % High-pass-filter data
    disp('HP-filtering');
    filtSweepSepData=fftFilter(filtSweepSepData',Fs,lower_cutoff,2);
    filtSweepSepData=filtSweepSepData';
    disp('Done HP-filtering');
end

LFPbySweep=filtSweepSepData;
