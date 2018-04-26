function [LFPbySweep,Fs,p,LFPspecgram]=makeWaveletSpecgram_wrapper(expt,fileInd,LFPbySweep)

Fs=32000;
downSampFactor=1;
noLEDcond=[5];
physCh=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];
dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
useTheseStimConds=1:9;
specgramFromThisChInd=4; % was 2

daqFileNames=expt.files.names(fileInd);
disp('Using these daq files');
disp(daqFileNames);

if isempty(LFPbySweep)
    % Get LFP
    [LFPbySweep,Fs,allCompleteSweeps]=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physCh(specgramFromThisChInd));
    if any(size(LFPbySweep{1})==0)
        disp('No data in this daq file.');
    end
    
    % Match acquired LFP sweeps to LED conditions
    ledConds=expt.sweeps.led(ismember(expt.sweeps.fileInd,fileInd));
    if length(ledConds)~=size(LFPbySweep{1},1)
        disp('ledConds and acquired LFP sweeps do not match');
        return
    end
    
    % Get stimconds
    stimConds=expt.sweeps.stimcond(ismember(expt.sweeps.fileInd,fileInd));
    if length(stimConds)~=size(LFPbySweep{1},1)
        disp('stimConds and acquired LFP sweeps do not match');
        return
    end
    
    % Get LFP sweeps without LED on
    for i=1:length(LFPbySweep)
        temp=LFPbySweep{i};
        %     LFPbySweep{i}=temp(ledConds==noLEDcond & ismember(stimConds,useTheseStimConds),:);
        %     LFPbySweep{i}=temp(ledConds==noLEDcond,:);
        LFPbySweep{i}=temp(ismember(ledConds,noLEDcond) & ismember(stimConds,useTheseStimConds),:);
    end
else
end

% Get average LFP
for i=1:length(LFPbySweep)
    avLFPbySweep{i}=mean(LFPbySweep{i},1);
    % Shift
    avLFPbySweep{i}=avLFPbySweep{i}-mean(avLFPbySweep{i});
end

% Display average LFP
figure();
cumOffset=0;
for i=1:length(avLFPbySweep)
    o=max(avLFPbySweep{i});
    cumOffset=cumOffset+o;
    plot(0:1/Fs:(length(avLFPbySweep{i})-1)*(1/Fs),avLFPbySweep{i}-cumOffset);
    hold on;
end

p=[];
LFPspecgram=[];
% [p,LFPspecgram]=makeWaveletSpecgram(LFPbySweep{1},1:size(LFPbySweep{1},1),Fs);