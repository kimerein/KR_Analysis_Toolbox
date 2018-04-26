function [f,avSpec,allSpecs]=makePowerSpectrum_wrapper(expt,fileInd,useTheseTrials)

timeWindow=[1.3 3];
Fs=32000;
downSampFactor=1;
useLEDcond=[1.050 2.050 4.050 6.050 8.050 10.050 12.050 14.050 16.050 18.050 20.050 30.050 40.050 50.050 60.050];
useTheseStimConds=1;
physChs=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];
choosePhysCh=2;
physCh=physChs(choosePhysCh);
% dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
dataDir='W:\New Acquisition Computer\';

daqFileNames=expt.files.names(fileInd);
disp('Using these daq files');
disp(daqFileNames);

% Get LFP
[LFPbySweep,Fs,allCompleteSweeps]=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physCh);
if any(size(LFPbySweep{1})==0)
    disp('No data in this daq file.');
end

% Match acquired LFP sweeps to LED conditions
ledConds=expt.sweeps.led(ismember(expt.sweeps.fileInd,fileInd));
if length(ledConds)~=size(LFPbySweep{1},1)
    disp('ledConds and acquired LFP sweeps do not match');
    dataInt=[];
    return
end

% Get stimconds
stimConds=expt.sweeps.stimcond(ismember(expt.sweeps.fileInd,fileInd));

% Get LFP sweeps with chosen LED cond
for i=1:length(LFPbySweep)
    temp=LFPbySweep{i};
    LFPbySweep{i}=temp(ismember(ledConds,useLEDcond) & ismember(stimConds,useTheseStimConds),:);
end

if strcmp(useTheseTrials,'all')
else
    % Filter according to useTheseTrials
    temp=LFPbySweep{1};
    if size(temp,1)~=length(useTheseTrials)
        disp('Error: useTheseTrials length does not match number of filtered LFP trials');
        dataInt=[];
        return
    end
    
    % Get useTheseTrials LFPSweeps
    for i=1:length(LFPbySweep)
        temp=LFPbySweep{i};
        LFPbySweep{i}=temp(useTheseTrials,:);
    end
end
    
% Display average LFP
figure(); 
xpoints=0:1/Fs:(size(LFPbySweep{1},2)-1)*(1/Fs);
plot(xpoints,mean(LFPbySweep{1},1));

% Plot power spectrum
temp=LFPbySweep{1};
[f,avSpec,allSpecs]=makePowerSpectrum(temp(:,xpoints>timeWindow(1) & xpoints<timeWindow(2)),Fs);


