function [Fs,ledDatabySweep,completeSweeps]=readTrialTriggeredLED(currFile,dataDir,daqFileName,Fs,downSampleFactor,ledChannel)
% Wrapper for daqread
% Reads the daq file containing the LED data in daqFileName
% 
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
data=daqread(strcat(dataDir,daqFileName),'Channels',[ledChannel]);
% physData=data(:,physChannel);
% photoData=data(:,photoChannel);
% ledData=data(:,ledChannel);
ledData=data;

% Parse ledData into sweeps
sweepSeps=find(isnan(ledData));
%disp('Loading DAQ LED data');
if numel(sweepSeps)~=0
    ledSepData=zeros(length(sweepSeps)+1,sweepSeps(1)-1);
    completeSweeps=zeros(length(sweepSeps)+1,1);
    ledSepData(1,:)=ledData(1:sweepSeps(1)-1)';
    disp(currFile);
    completeSweeps(1)=1;
    for i=1:length(sweepSeps)
        if mod(i,100)==0
            %disp('Loading DAQ LED data');
        end
        if i==length(sweepSeps)
            if length(ledData)-sweepSeps(end)==length(ledSepData(1,:))
                ledSepData(i+1,:)=ledData(sweepSeps(i)+1:end)';
                disp(currFile);
                completeSweeps(i+1)=1;
            else
                ledSepData(i+1,:)=[ledData(sweepSeps(i)+1:end); zeros(length(ledSepData(1,:))-(length(ledData)-sweepSeps(end)),1)]';
                disp(currFile);
                completeSweeps(i+1)=0;
            end
        else
            ledSepData(i+1,:)=ledData(sweepSeps(i)+1:sweepSeps(i+1)-1)';
            disp(currFile);
            completeSweeps(i+1)=1;
        end
    end
else
    completeSweeps=1;
    ledSepData=ledData';
    disp(currFile);
end

if downSampleFactor~=1
    currFs=Fs;
    [ledDatabySweep lFs]=downSamp(ledSepData,currFs,downSampleFactor);
else
    ledDatabySweep=ledSepData;
end